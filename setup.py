# -*- coding:utf-8 -*-
from distutils.core import setup, Extension
import os

# Compile the list of packages available, because distutils doesn't have
# an easy way to do this.
packages, data_files = [], []
root_dir = os.path.dirname(__file__)
if root_dir:
    os.chdir(root_dir)

'''注册Extension, 通过运行 python setup.py build 来编译模块'''
geofix_all = Extension('geolocationfix.geofix_all',
                    sources = ['geofix_all.c'])

setup(name='geolocationfix',
      description='An extensible geolocation fix tool',
      author='wfxiang08',
      version = '1.1',
      author_email='wfxiang08@gmail.com',
      url='http://www.bitbucket.org/ubernostrum/django-registration/wiki/',
      download_url='http://www.bitbucket.org/ubernostrum/django-registration/get/v0.7.gz',
      packages=['geolocationfix'],
      package_data={'geolocationfix': ['Gps2Sogou.txt', 'Sogou2Wgs.txt']},
      ext_modules = [geofix_all]
      )
