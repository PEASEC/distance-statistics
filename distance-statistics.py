import json
import os
import pickle
from typing import TextIO

import geopy.distance
import shapely.geometry
import geopandas as gp
from pandas import DataFrame
from shapely.geometry import shape
import pandas as pd
import numpy as np
import random
from OSMPythonTools.overpass import Overpass, OverpassResult
from OSMPythonTools.cachingStrategy import CachingStrategy, JSON
from os.path import exists
from geopy import distance

CachingStrategy.use(JSON, cacheDir='.osm_cache')

overpass_endpoint = "http://overpass-api.de/api/"
#overpass_endpoint = "http://127.0.0.1:12345/api/"


# Helper filepath functions #

def unprocessed_farmyards_filepath(state: str, county: str):
    return f"tmp/unprocessed_farmyards/{state}-{county}-unprocessed_farmyards.json"


def merged_farmyards_filepath(state: str, county: str):
    return f"tmp/merged_farmyards/{state}-{county}-merged_farmyards.json"


def representative_buildings_filepath(state: str, county: str, idx: int):
    return f"tmp/representative_buildings/{state}-{county}-{idx}-representative_buildings.json"


def buildings_filepath(state: str, county: str):
    return f"tmp/buildings/{state}-{county}-buildings.json"


def distance_matrix_filepath(state: str, county: str):
    return f"distance_matrices/{state}-{county}-distance_matrix.csv"


def min_dist_filepath(state: str, county: str):
    return f"evaluation/minimum_distances/{state}-{county}-minimum_distances.csv"


def lowest_five_dists_filepath(state: str, county: str):
    return f"evaluation/minimum_five_dists/{state}-{county}-minimum_five_dists.csv"


def count_of_farms_in_radius_filepath(state: str, county: str):
    return f"evaluation/count_of_farms_in_radius/{state}-{county}-count_of_farms_in_radius.csv"


# Main functions separated into numbered steps #

def step0_get_counties(state: str):
    overpass = Overpass(endpoint=overpass_endpoint)

    timeout = 1000
    out = 'json'

    try:
        # Get sub-states on admin-level 6
        query = f"""
area[admin_level="4"][boundary="administrative"][name="{state}"]->.searchArea;
(
  rel(area.searchArea)[boundary="administrative"]["admin_level"="6"];
);
map_to_area;
out geom;"""
    except:
        print('Exception during overpass query to get county.')
        return None

    counties_result = overpass.query(query, timeout=timeout, out=out)
    return counties_result


def step1_get_farm_areas(state: str):
    counties_result = step0_get_counties(state)

    if counties_result.areas() is None:
        return

    for i in range(len(counties_result.areas())):
        idx = counties_result.areas()[i].id()
        county = counties_result.areas()[i].tag('name')
        county = county.replace("/", ":")

        filepath = unprocessed_farmyards_filepath(state, county)
        if os.path.exists(filepath):
            print(f"skip: {filepath} exists already")
            continue

        timeout = 1000
        out = 'json'
        query = \
            f"""
         way(area:{idx})["landuse"="farmyard"]->.all_farmyards;
         .all_farmyards out geom;
        """
        try:
            overpass = Overpass(endpoint=overpass_endpoint)
            result = overpass.query(query, timeout=timeout, out=out)

            if result.ways() is None:
                print(f"skip - {state}-{county} has no ways (areas) tagged with farmyard!")
                continue

            shapes = []
            for area in result.ways():
                shapes.append(area.geometry())

            gs = gp.GeoSeries(data=shapes, crs="EPSG:4326")
            gs.to_file(filepath, driver="GeoJSON")
        except:
            print(f"Exception during overpass query to get farmyard of county {county}")


def step2_filter_neighbored_areas(state: str):
    print(f"  Filter neighbored areas for state '{state}'")
    counties_result = step0_get_counties(state)
    if counties_result.countWays() is None:
        return None

    for county_i in range(len(counties_result.areas())):
        county = counties_result.areas()[county_i].tag('name')
        county = county.replace("/", ":")
        filepath_merged_farmyards = merged_farmyards_filepath(state, county)
        if os.path.exists(filepath_merged_farmyards):
            print(f"    skip {state}-{county} - exists already")
            continue

        filepath_unprocessed_farmyards = unprocessed_farmyards_filepath(state, county)
        if not os.path.exists(filepath_unprocessed_farmyards):
            print(f"    skip - {filepath_unprocessed_farmyards} does not exist")
            continue

        print(f"    County '{county}'")
        gs = gp.GeoSeries.from_file(filepath_unprocessed_farmyards)

        # Check distances between
        result = []
        filtered_id = []  # keep track of already sorted out areas
        filter_distance = 0.30  # filter distance in kilometers

        # iterate over all shapes and filter out areas based on their distance next to each other
        new_polys = gs.array.tolist()

        filter_areas = []
        for i in range(len(new_polys)):
            if not hasattr(new_polys[i], "exterior"):
                print(f"      filter {i}")
                filter_areas.append(i)
        filtered = [ele for idx, ele in enumerate(new_polys) if idx not in filter_areas]
        new_polys = filtered

        for i in range(len(new_polys)):
            if i in filtered_id:
                continue
            poly = new_polys[i]
            # check distance with all remaining polys
            poly_center = (poly.centroid.y, poly.centroid.x)
            if i % 100 == 0:
                print(f"      i = {i}")
            for j in range(len(new_polys)):
                if j in filtered_id or i == j:
                    continue
                poly_test = new_polys[j]
                poly_test_center = (poly_test.centroid.y, poly_test.centroid.x)
                if distance.distance(poly_center, poly_test_center).km < filter_distance:
                    # merge the polygons
                    filtered_id.append(j)
                    points = shapely.geometry.MultiPoint(new_polys[i].exterior.coords[:-1]
                                                         + new_polys[j].exterior.coords[:-1])

                    new_polys[i] = points.convex_hull
                    poly = new_polys[i]  # if next iteration tries to access poly ...

        #  remove filtered_ids from gs
        result = [ele for idx, ele in enumerate(new_polys) if idx not in filtered_id]

        gs = gp.GeoSeries(data=result, crs="EPSG:4326")  # evtl. hier noch x-y vertauschen
        gs.to_file(filepath_merged_farmyards, driver="GeoJSON")

        print(f"      len of result is {len(result)} - original {len(new_polys)}")


def step3_get_representative_buildings(state: str):
    print(f"  Get representative buildings for state '{state}'")
    counties_result = step0_get_counties(state)
    if counties_result.countWays() is None:
        return None

    for county_i in range(len(counties_result.areas())):
        county = counties_result.areas()[county_i].tag('name')
        county = county.replace("/", ":")

        filepath_buildings = buildings_filepath(state, county)
        if os.path.exists(filepath_buildings):
            print(f"    skip - {filepath_buildings} exists already")
            continue

        filepath_merged_farmyards = merged_farmyards_filepath(state, county)
        if not os.path.exists(filepath_merged_farmyards):
            print(f"    skip - {filepath_merged_farmyards} does not exist")
            continue

        gs = gp.GeoSeries.from_file(filepath_merged_farmyards, driver="GeoJSON")
        if gs.size == 0:
            continue

        print(f"    get all representative buildings for county '{county}'")
        polys = gs.array.tolist()  # polygons of county

        county_buildings = []

        # request overpass for all buildings in a polygon area:
        i: int
        for i, poly in enumerate(polys):

            poly_str = ""
            if not (len(poly.exterior.coords) % 2 == 0):
                poly.exterior.coords = poly.exterior.coords[:-1]
            for coord in poly.exterior.coords:
                poly_str += str(coord[1]) + " " + str(coord[0]) + "  "
            poly_str = poly_str.rstrip()

            print(poly_str)
            timeout = 1000
            out = 'json'
            query = \
                f"""
                way(poly:"{poly_str}")[building];
                out geom;
            """
            try:
                overpass = Overpass(endpoint=overpass_endpoint)
                result = overpass.query(query, timeout=timeout, out=out)
            except:
                print(f"    Exception during overpass query to get farmyard of county {county}")
                continue

            if result.ways() is None:
                print(f"    skip - {state}-{county}-{i} has no building at all")
                continue

            shapes = []
            for area in result.ways():
                shapes.append(area.geometry())

            # store intermediate result
            gs = gp.GeoSeries(data=shapes, crs="EPSG:4326")
            filepath_representative_buildings = representative_buildings_filepath(state, county, i)
            gs.to_file(filepath_representative_buildings, driver="GeoJSON")

            # save one random building
            random_shape_idx = random.randrange(len(shapes))
            county_buildings.append(shapes[random_shape_idx])

        if len(county_buildings) == 0:
            print(f"    skip - {state}-{county} has no buildings on farmyard")
            continue

        gs = gp.GeoSeries(data=county_buildings, crs="EPSG:4326")
        gs.to_file(filepath_buildings, driver="GeoJSON")


def step4_calculate_distance_matrix(state: str):
    print(f"  Calculate distance matrices for state '{state}'")
    counties_result = step0_get_counties(state)
    if counties_result.countWays() is None:
        return None

    for county_i in range(len(counties_result.areas())):
        county = counties_result.areas()[county_i].tag('name')
        county = county.replace("/", ":")

        filepath_buildings = buildings_filepath(state, county)
        if not os.path.exists(filepath_buildings):
            print(f"    Could not load file {filepath_buildings}")
            continue

        gs = gp.GeoSeries.from_file(filepath_buildings, driver="GeoJSON")
        polys = gs.array.tolist()  # polygons of buildings

        all_points = []
        for poly in polys:
            all_points.append((poly.centroid.y, poly.centroid.x))  # for geopy: use (lat, long), i.e., (y, x)

        # use pickle, as calculation can be a bit tense in terms of time and memory
        # ... so the calculation can proceed at a previous state
        filepath_dump_matrix = f"tmp/pickle/.dist_mat_{state}_{county}.dump"
        filepath_dump_idx = f"tmp/pickle/.dist_mat_{state}_{county}_i.dump"
        if os.path.exists(filepath_dump_idx) and os.path.exists(filepath_dump_matrix):
            with open(filepath_dump_idx, "rb") as pickle_off:
                dump_i = int(pickle_off.readline(100).strip())
            with open(filepath_dump_matrix, "rb") as pickle_off:
                dist_matrix = pickle.load(pickle_off)
        else:
            dump_i = 0
            dist_matrix = np.empty(shape=(len(all_points), len(all_points)), dtype=np.float16)

        for i in range(dump_i, len(all_points)):
            point_i = all_points[i]
            for j in range(i, len(all_points)):
                dist_matrix[i, j] = distance.distance(point_i, all_points[j]).km
            if (i % 100) == 0:
                print(f"      dump after row {i}/{len(all_points)}")
                dist_matrix.dump(filepath_dump_matrix)
                with open(filepath_dump_idx, "w") as pickle_i_file:
                    pickle_i_file.write(str(i))

        for i in range(0, len(all_points)):
            for j in range(0, i):
                dist_matrix[i, j] = dist_matrix[j, i]

        filepath_matrix = distance_matrix_filepath(state, county)
        gdf = gp.GeoDataFrame(dist_matrix, index=all_points, columns=all_points)
        gdf.to_csv(filepath_matrix)


def step5_create_evaluation_files(state: str):
    print(f"  Evaluate distance matrices of {state}")

    counties_result = step0_get_counties(state)
    if counties_result.countWays() is None:
        print(f"    error - could not find valid counties for '{state}'")
        return None

    for county_i in range(len(counties_result.areas())):
        county = counties_result.areas()[county_i].tag('name')
        county = county.replace("/", ":")

        filepath_distance_matrix = distance_matrix_filepath(state, county)
        if not os.path.exists(filepath_distance_matrix):
            print(f"    error - could not load '{filepath_distance_matrix}'")
            continue

        df = pd.read_csv(filepath_distance_matrix, index_col=0)
        if len(df.columns) == 1:
            print(f"    has just one entry")
            continue

        minimum_distances_count = 5

        minimum_distances_count_dict: dict = {}
        filepath_min_dist = min_dist_filepath(state, county)
        filepath_lowest_five = lowest_five_dists_filepath(state, county)
        if not exists(filepath_min_dist) or not exists(filepath_lowest_five):
            for point, distances in df.iteritems():
                minimum_distances_count_dict[point] = distances.nsmallest(minimum_distances_count + 1)[1:].array

        # part 1: minimum distance stats
        if not os.path.exists(filepath_min_dist):
            minimum_distances = []
            for point in minimum_distances_count_dict:
                minimum_distances.append(minimum_distances_count_dict[point][0])

            minimum_distances = pd.DataFrame(minimum_distances)
            minimum_distances.to_csv(filepath_min_dist)

        # part 2: minimum 5 distances (kind of superset of part 1)
        if not exists(filepath_lowest_five):
            minimum_five_dists = []
            for point in minimum_distances_count_dict:
                minimum_five_dists.append(minimum_distances_count_dict[point])

            minimum_five_dists = pd.DataFrame(minimum_five_dists)
            minimum_five_dists.to_csv(filepath_lowest_five)

        # part 3: count of farms in radius of x km | x e {1,3,5,7,9} km
        filepath_farms_in_radius = count_of_farms_in_radius_filepath(state, county)
        if not exists(filepath_farms_in_radius):
            count_of_farms_in_radius = []
            radius_list = [1, 2, 3, 4, 5, 7, 9]
            for point, distances in df.iteritems():
                radius_neighbours_list = []
                dist_array = distances.to_numpy()
                for radius in radius_list:
                    radius_neighbour_counter = np.count_nonzero(dist_array < radius) - 1  # we won't count itself (0)
                    radius_neighbours_list.append(radius_neighbour_counter)
                count_of_farms_in_radius.append(radius_neighbours_list)

            count_of_farms_in_radius = pd.DataFrame(count_of_farms_in_radius, columns=radius_list)
            count_of_farms_in_radius.to_csv(filepath_farms_in_radius)


def step6_write_stats_summary(state: str, file: TextIO):
    print(f"  Generate stats for {state}")

    file.write(f"\n# {state}\n")

    counties_result = step0_get_counties(state)
    if counties_result.countWays() is None:
        print(f"    error - could not find valid counties for '{state}'")
        return None

    complete_data = []

    for county_i in range(len(counties_result.areas())):
        county_original = counties_result.areas()[county_i].tag('name')
        county = county_original.replace("/", ":")

        county_data = {'state': state, 'county': county}

        filepath = distance_matrix_filepath(state, county)
        if not os.path.exists(filepath):
            print(f"    error - could not load '{filepath}'")
            continue

        df = pd.read_csv(filepath, index_col=0)
        if len(df.columns) == 1:
            print(f"    '{county}' has just one entry")
            continue

        file.write(f"\n## {county}\n\n")

        filepath_min_dist = min_dist_filepath(state, county)
        filepath_lowest_five = lowest_five_dists_filepath(state, county)
        filepath_count_in_radius = count_of_farms_in_radius_filepath(state, county)
        min_dists = pd.read_csv(filepath_min_dist, index_col=0)
        file.write(f"    count of farms = {len(min_dists)}  \n\n")

        file.write(f"    min.min  = {min_dists.min().to_string(index=False)}  \n"
                   f"    min.max  = {min_dists.max().to_string(index=False)}  \n"
                   f"    min.mean = {min_dists.mean().to_string(index=False)}  \n"
                   f"    min.std  = {min_dists.std().to_string(index=False)}  \n\n")

        lowest_five = pd.read_csv(filepath_lowest_five, index_col=0)
        for i in lowest_five.columns[1::]:
            file.write(f"    lowest_five[{i}].min  = {lowest_five[i].min()}  \n"
                       f"    lowest_five[{i}].max  = {lowest_five[i].max()}  \n"
                       f"    lowest_five[{i}].mean = {lowest_five[i].mean()}  \n"
                       f"    lowest_five[{i}].std  = {lowest_five[i].std()}  \n\n")
            county_data[f'lowest_five_{i}_min'] = lowest_five[i].min()
            county_data[f'lowest_five_{i}_max'] = lowest_five[i].max()
            county_data[f'lowest_five_{i}_mean'] = lowest_five[i].mean()
            county_data[f'lowest_five_{i}_std'] = lowest_five[i].std()

        count_in_radius = pd.read_csv(filepath_count_in_radius, index_col=0)
        for i in count_in_radius.columns:
            file.write(f"    count_in_radius[{i}km].min  = {count_in_radius[i].min()}  \n"
                       f"    count_in_radius[{i}km].max  = {count_in_radius[i].max()}  \n"
                       f"    count_in_radius[{i}km].mean = {count_in_radius[i].mean()}  \n"
                       f"    count_in_radius[{i}km].std  = {count_in_radius[i].std()}  \n")
            county_data[f'count_in_radius_{i}km_min'] = count_in_radius[i].min()
            county_data[f'count_in_radius_{i}km_max'] = count_in_radius[i].max()
            county_data[f'count_in_radius_{i}km_mean'] = count_in_radius[i].mean()
            county_data[f'count_in_radius_{i}km_std'] = count_in_radius[i].std()

        overpass = Overpass(endpoint=overpass_endpoint)
        query = f"""
        area["admin_level"=4]["name"="{state}"]->.a;
        rel[name="{county_original}"][type=boundary]["admin_level"=6](area.a);
        out geom;"""

        timeout = 1000
        out = 'json'
        result = overpass.query(query, timeout=timeout, out=out)

        geometry = result.relations()[0].geometry()
        county_data['geometry'] = geometry

        complete_data.append(county_data)

    gdf = gp.GeoDataFrame(data=complete_data, crs="EPSG:4326")
    gdf.to_file(f"evaluation/complete_statistics/{state}-complete_statistics.geojson", driver="GeoJSON")


def main():
    states = [
        "Baden-Württemberg",
        "Bayern",
        "Berlin ",
        "Brandenburg",
        "Bremen ",
        "Hamburg ",
        "Hessen",
        "Mecklenburg-Vorpommern",
        "Niedersachsen",
        "Nordrhein-Westfalen",
        "Rheinland-Pfalz",
        "Saarland",
        "Sachsen",
        "Sachsen-Anhalt",
        "Schleswig-Holstein",
        "Thüringen"
    ]

    for state in states:
        print(f"# {state}\n")
        step1_get_farm_areas(state)
        step2_filter_neighbored_areas(state)
        step3_get_representative_buildings(state)
        step4_calculate_distance_matrix(state)
        step5_create_evaluation_files(state)

    filepath_stats = "stats_summary.md"
    overwrite = True
    if not overwrite and exists(filepath_stats):
        print(f"File {filepath_stats} exists already. Does not generate new stats file.")
    else:
        with open(filepath_stats, "a") as file:
            for state in states:
                step6_write_stats_summary(state, file)


if __name__ == '__main__':
    main()
