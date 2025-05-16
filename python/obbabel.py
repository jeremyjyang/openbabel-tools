#!/usr/bin/env python
###
# http://openbabel.org/docs/FileFormats/Overview.html
# Also list all formats with "obabel -L formats"
# http://openbabel.org/docs/Command-line_tools/babel.html
# http://openbabel.org/docs/UseTheLibrary/Python_Pybel.html
###
import sys,os,argparse,logging
import pandas as pd

import openbabel
from openbabel import pybel

#############################################################################
def Convert(ifile, ofile, iformat, oformat, nmax):
  mol = openbabel.OBMol(); n_out=0;
  obconv = openbabel.OBConversion()

  if not iformat: iformat = obconv.FormatFromExt(ifile)
  if not oformat: oformat = obconv.FormatFromExt(ofile)

  if not iformat or not obconv.SetInFormat(iformat):
    logging.error(f"Cannot read input format: {iformat}")
  if not oformat or not obconv.SetOutFormat(oformat):
    logging.error(f"Cannot write output format: {oformat}")

  while True:
    mol.Clear()
    if n_out==0:
      ok = obconv.ReadFile(mol, ifile)
      ifs = obconv.GetInStream()
    else:
      ok = obconv.Read(mol, ifs)
    if not ok or mol.Empty():
      break
    if nmax and n_out==nmax:
      break
    if n_out==0:
      ok = obconv.WriteFile(mol, ofile)
      ofs = obconv.GetOutStream()
    else:
      ok = obconv.Write(mol, ofs)
    if not ok:
      logging.warning("Warning: write problem.")
    n_out+=1
  logging.info(f"mols processed: {n_out}")

#############################################################################
def ListInputFormats():
  df = pd.DataFrame({'ext':pybel.informats.keys(), 'description':pybel.informats.values()})
  df.to_csv(sys.stdout, sep="\t", index=False)
  logging.info(f"Input formats: {df.shape[0]}")

#############################################################################
def ListOutputFormats():
  df = pd.DataFrame({'ext':pybel.outformats.keys(), 'description':pybel.outformats.values()})
  df.to_csv(sys.stdout, sep="\t", index=False)
  logging.info(f"Output formats: {df.shape[0]}")

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='OpenBabel chemical file converter')
  ops = [ "list_input_formats", "list_output_formats", "convert", ]
  parser.add_argument("op", choices=ops, help='OPERATION (select one)')
  parser.add_argument("--i", dest="ifile", help="input file")
  parser.add_argument("--o", dest="ofile", help="output (TSV)")
  parser.add_argument("--iformat")
  parser.add_argument("--oformat")
  parser.add_argument("--nmax", type=int, default=None)
  parser.add_argument("-v","--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  if args.op == "list_input_formats":
    ListInputFormats()
 
  elif args.op == "list_output_formats":
    ListOutputFormats()
 
  elif args.op == "convert":
    Convert(args.ifile, args.ofile, args.iformat, args.oformat, args.nmax)

  else:
    parser.error(f'Invalid operation: {args.op}')
