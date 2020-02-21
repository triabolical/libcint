import sys
import csv
import numpy
import cProfile
import ctypes
import warnings
import pstats
import math
import time
import input
import numpy as np
import os.path
from pyscf.df import addons
import scipy
from scipy.linalg import fractional_matrix_power
from pyscf import gto, scf, lib, ao2mo, __config__, df, pbc

atom = [];
atom.append(['H', (0.5, 0.5, 0.5)])
atom.append(['H', (0.5, 0.5, 0.5)])
mol = gto.mole.Mole()
mol.atom = atom
mol.build()
M = mol.intor('int1e_nuc_sph')
print(M)
