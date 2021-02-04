EO2RIVER / EO2FM

* working python tool (command-line) to generate river network with cross-sections given polygon (servers-side Python + Google Earth Engine)

* working Python library (API, Python Notebook)

* generate cross-sections using fusion of different datasets (from discharge/width/slope, water occurrence, elevation model, ICESat-2 tracks)
    
    What discharge source(s) to use (GLOFFIS, HydroRIVERS, GRDC, etc.)?

* improve algorithm (branch pruning, select main river)

* generate high-resolution water occurrence (better than JRC)

* make use of OpenStreetMap data

* make the tool work for larger areas (parallel processing)

* prepare a number of test demo cases showing functionality of the tool

    * Brahmaputra (braided patters every year)
    * Niger (arid, large water level changes)
    * Amazon tributaries (meandering)
    * Rhine (cloud challengind, no intra-annual changes)
    * Small River (example Phom Penh, mix OSM, JRC, HydroRivers, etc.)

* integration with HydroLib project, dhydamo - QGIS plugin

* integration with Blue Earth (which parts?)

* click point - select upstream rivers (taking catchments into account) - ???
