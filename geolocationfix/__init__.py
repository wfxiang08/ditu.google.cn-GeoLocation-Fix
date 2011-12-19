# -*- coding:utf-8 -*-
import os
_modulePath = os.path.dirname(__file__)
destFileName = os.path.join(_modulePath, "geofix_all.so")

if not os.path.exists(destFileName):
    raise Exception('Error', "Configuration File Gps2Sogou.tar.gz not found")

_filePath1 = os.path.abspath(os.path.join(_modulePath, "Sogou2Wgs.txt"))
_filePath2 = os.path.abspath(os.path.join(_modulePath, "Gps2Sogou.txt"))
if not os.path.exists(_filePath1):
    if not os.path.exists(os.path.join(_modulePath, "Sogou2Wgs.tar.gz")):
        raise Exception('Error', "Configuration File Sogou2Wgs.tar.gz not found")

if not os.path.exists(_filePath2):
    if not os.path.exists(os.path.join(_modulePath, "Gps2Sogou.tar.gz")):
        raise Exception('Error', "Configuration File Gps2Sogou.tar.gz not found")
