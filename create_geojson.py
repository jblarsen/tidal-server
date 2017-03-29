#!/bin/env python
# -*- coding: utf-8 -*-
"""
Convert data from DMI Greenland tidal station report to GeoJSON.
"""
import glob
import sys
import json
import os.path
import xml.etree.cElementTree as ET

def read_station(filename):
    """Read in DMI station names and positions."""
    data = []
    with open(filename, 'r') as ifile:
        sid = os.path.basename(filename)
        sid = sid.split('.')[0]
        while True:
            try:
                line = ifile.next()
                if line.startswith('Station / Station:'):
                    name = line.split(':')[1].strip()
                if line.startswith('Bredde / Latitude:'):
                    lat = line.split(':')[1]
                    lat = float(lat.split('(')[1].split(')')[0])
                if line.startswith('Længde / Longitude:'):
                    lon = line.split(':')[1]
                    lon = float(lon.split('(')[1].split(')')[0])
            except StopIteration:
                break
    return (sid, name, lat, lon)

def main():
    indir = sys.argv[1]

    stations = []
    for infile in glob.glob(indir + '/*.tmp.txt'):
        station = read_station(infile)
        stations.append(station)

    feature_coll = { "type": "FeatureCollection",
                     "features": []}

    for station in stations:
        sid, name, lat, lon = station
        feature = {"type": "Feature",
                   "geometry": {"type": "Point", "coordinates": [lon, lat]},
                   "properties": {"id": sid,
                                  "name": name}
                   }
        feature_coll["features"].append(feature)
  
    print json.dumps(feature_coll, indent=2)

if __name__ == '__main__':
    main()
