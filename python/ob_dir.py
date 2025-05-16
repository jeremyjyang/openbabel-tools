#!/usr/bin/env python

import openbabel

tokens = dir(openbabel)

for token in tokens:
  print('openbabel.%s'%token)
