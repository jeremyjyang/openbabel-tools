#!/usr/bin/env python
#############################################################################
import os,sys
import openbabel


#############################################################################
def OBCansmi(smiles, isomeric=False):
  obc = openbabel.OBConversion()
  mol = openbabel.OBMol()
  obc.SetInFormat("smi")
  obc.ReadString(mol, smiles)
  if not isomeric:
    obc.AddOption("i", obc.OUTOPTIONS)
  obc.SetOutFormat("can")
  cansmi = obc.WriteString(mol, 1)
  return cansmi

#############################################################################
