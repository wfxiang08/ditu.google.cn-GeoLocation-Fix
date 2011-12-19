# -*- coding:utf-8 -*-
import geofix_all #IGNORE:100
import logging
import os
import re
import time

class GeoFixAll(object):
    
    PATTERN = re.compile("\\d+", re.DOTALL)    
    def __init__(self):
        _modulePath = os.path.dirname(__file__)
        _filePath1 = os.path.abspath(os.path.join(_modulePath, "Sogou2Wgs.txt"))
        _filePath2 = os.path.abspath(os.path.join(_modulePath, "Gps2Sogou.txt"))
        geofix_all.loadAllData(_filePath1, _filePath2)

    def sogouToGoogleMars(self, x, y):
        '''
            @param x, y: x coordinate in Sogou Map System(和longitude差不多)
            返回tuple(longitude, latitude)
        '''        
        result = geofix_all.sogouToGoogleMars(x, y)
        return result

    def gps2Sogou(self, lng, lat):
        '''
            @param lng, lat: gps in degree 
            @return: (x, y)
        '''        
        result = geofix_all.gpsToSogou(lng, lat)
        return result

    def gpsToGoogleMars(self, lng, lat):
        result = geofix_all.gpsToGoogleMars(lng, lat)
        return result
                 
'''创建静态变量'''
_geoFixAll = GeoFixAll()

def sogouToGoogleMars(x, y):
    '''
        @param x, y: x coordinate in Sogou Map System(和longitude差不多)
        返回tuple(longitude, latitude)
    '''
    return _geoFixAll.sogouToGoogleMars(x, y)

def gpsToGoogleMars(x, y):
    '''
        @param lng, lat: gps in degree 
        @return: (lng, lat) in Google Map system
    '''    
    return _geoFixAll.gpsToGoogleMars(x, y)


def testGpsToMars(y, x, expY, expX):
    result = gpsToGoogleMars(x, y)
    print "(%10.5f, %10.5f) ==> (%10.5f, %10.5f) vs. (%10.5f, %10.5f)" % (y, x, result[1], result[0], expY, expX)

if __name__ == "__main__":
    result = sogouToGoogleMars(11349877.565125303, 4262240.989378193)
    print "Act: (%f, %f)  <===> Exp: 101.95832890731322,35.904431065652595" % result
    testGpsToMars(40.00688, 116.34395, 0, 0)