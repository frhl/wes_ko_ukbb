#!/usr/bin/env python3

import re

def extract_from_header(path, regexes = ['ID','age','PC+','ukbb', 'sex','array'], delim = '\t'):
    '''Extract phenotypes from the header by removing regular expression matches (e.g. covaraties)'''
    combined = "(" + ")|(".join(regexes) + ")"
    infile = open(path,'r')
    line = infile.readline().strip('\n').split(delim)
    line_clean = [l for l in line if not re.match(combined, l)]
    return line_clean


