#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os.path
import os
import sys
import pandas as pd

dataset_name = sys.argv[1]
timeout = int(sys.argv[2])

df = pd.read_csv(dataset_name, index_col=False)
df.columns = ["Time"]
column = df["Time"]
print df.head()



print "Mean: " + str(column.mean())
print "Median: " + str(column.median())
print "Timeout: " + str(column[column >=timeout].count())