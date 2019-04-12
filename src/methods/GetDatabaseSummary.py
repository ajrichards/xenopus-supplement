#!/usr/bin/python

from htsint.database import print_db_summary
printstr = print_db_summary()

logfile = open('dbsummary.log','w')
logfile.write(printstr)
logfile.close()
