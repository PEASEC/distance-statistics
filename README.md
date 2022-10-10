# Distance Statistics

## Purpose and Motivation

This project consists of scripts to calculate and evaluate distances between specific building types, e.g., farms
It was created to analyze the geographic distribution of farms in Germany, but could be adapted to calculate other 
objects in arbitrary regions.

## Dependencies

* python3
* python3 packages (installation inside venv, see 'Project Setup')
  * OSMPythonTools (https://github.com/mocnik-science/osm-python-tools)
  * Shapely
  * Pandas
  * Geopandas
  * Numpy
  * Geopy
  * folium
  * Mapclassify

## Project Setup

Create local venv:

    $ python3 -m venv ./venv

Activate venv:

    $ source ./venv/bin/activate

Install dependencies:

    $ pip install OSMPythonTools shapely pandas numpy geopy folium mapclassify

## Usage 

### 1. Calculate distance statistics

The `distance-statistics.py` retrieves coordinates of farmyards in Germany via 
OverPass, filters the retrieved areas (remove farmyards without buildings), and
selecting one building as a representative farm. 
With this farm it calculates 
distance statistic files.
Retrieval and calculation is restricted to county-level of Germany. 

#### Configure

You can change parts of the code to select regions 

#### Execute

Activate venv:

    $ source ./venv/bin/activate

Run main script (*the script quite take some time*):

    $ python3 distance-statistics.py


This will create and filling following folders with results:  
* `distance_matrices/`
  CSV-files (county-level) containing N² farm-to-farm distance matrices
* `evaluation/`
  * `evaluation/complete_statistics/` 
    GeoJSON-files (state-level) containing concatenated statistics of 
    counties, including geometries of all farms
  * `evaluation/count_of_farms_in_radius/`
    CSV-files (county-level) containing count of neighbors in [1,2..9] km
    radius for each farm
  * `evaluation/minimum_distances/`
    CSV-files (county-level) containing the single minimum distance to another 
    farm (deprecated)
  * `evaluation/minimum_five_dists/`
    CSV-files (county-level) containing the five minimum distances to other
    farms (this supersedes minimum_distances)

And following directories for debugging and caching purpose (speed up on next run) 
* `tmp/`
  intermediate calculation artefacts, like retrieved farm-building geometries
* `.osm_cache/` 
  cached data from Overpass

#### Use own Overpass server

The public overpass server has rate limits.
For setting up an own OverPass-Server, see: 
https://wiki.openstreetmap.org/wiki/Overpass_API/Installation or 
https://hub.docker.com/r/wiktorn/overpass-api/

After setting up you have to change the following line at the begin of 
`distance-statistics.py` to the desired server address:

    overpass_endpoint = "http://overpass-api.de/api/"

For example to:

    overpass_endpoint = "http://127.0.0.1:12345/api/"

### 2. Evaluation

There exists two Jupyter Notebook files that support evaluation:
* `statistics.ipynb`
  Presents statistics, foremost count of farms with n [n=1,2..5] minimum 
  neighbours in m [m=1,2..5km] radius
* `clustering.ipynb`
  Creates clusters including plots

## Other tools / Debugging

### osmtogeojson: Converting cached overpass query data to GeoJson

See https://github.com/tyrasd/osmtogeojson for installation.

Create dir for intermediate files:

    mkdir json

Show existing cached data:

    ls .osm_cache

Use for example "overpass-50d9235b285ba1adddb68f96b061de39efa68db8": 

    export OVERPASS_CACHE=overpass-50d9235b285ba1adddb68f96b061de39efa68db8

Create formatting required by other tools like `osmtogeojson` with `jq`: 

    cat .osm_cache/${OVERPASS_CACHE} | jq > json/${OVERPASS_CACHE}.json

Remove top-level elements, as `osmtogeojson` cannot interpret them: 

    vim json/${OVERPASS_CACHE}.json

Convert Json to GeoJson: 

    osmtogeojson json/${OVERPASS_CACHE}.json > json/${OVERPASS_CACHE}.geojson

## Acknowledgments

This work was created at Science and Technology for Peace and Security (PEASEC), Technical University of Darmstadt, www.peasec.de, and supported by funds of the German Government’s Special Purpose Fund held at Landwirtschaftliche Rentenbank in the projects Geobox-II and AgriRegio.
* Contributors under those funds:
  * Franz Kuntke

## License
Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at your option.

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in `distance-statistics` by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any additional terms or conditions.
