#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Process Quantum Espresso pw.x output file into a CIF file. 
# 
# Ported from perl script qe2cif.pl by Alexandr Fonari
# Updated to read atomic positions from current version of QE pw.x output file. 
#
# Author: mikas.remeika@bk.tsukuba.ac.jp
# 
# Tested with 'vc-relax', QE v6.2.
# 
# Input parameters:
# 
# Which atomic positions to process: -sf (start and final) -int (intermediate).
# Intermediate positions do not include start and final positions.
# All positions will be written into one CIF file, labeles accordingly. 
#
# pw.x output file (captuted output to stdout). 

import argparse 
import re 
import math
import cmath
import numpy as np

def print_cif_header(a, b, c, alpha, beta, gamma, data_title):
    f_out.write("\n\ndata_{}\n".format(data_title))
    f_out.write("\n")
    f_out.write("_cell_length_a     {:16.14f}\n".format(a))
    f_out.write("_cell_length_b     {:16.14f}\n".format(b))
    f_out.write("_cell_length_c     {:16.14f}\n".format(c))
    f_out.write("_cell_angle_alpha  {:16.14f}\n".format(alpha))
    f_out.write("_cell_angle_beta   {:16.14f}\n".format(beta))
    f_out.write("_cell_angle_gamma  {:16.14f}\n".format(gamma))
    f_out.write("\n")
    f_out.write("loop_\n")
    f_out.write("_space_group_symop_id\n")
    f_out.write("_space_group_symop_operation_xyz\n")
    f_out.write("1 'x, y, z'\n")
    f_out.write("\n")
    f_out.write("loop_\n")
#    f_out.write("_atom_site_label\n")
    f_out.write("_atom_site_type_symbol\n")
    f_out.write("_atom_site_fract_x\n")
    f_out.write("_atom_site_fract_y\n")
    f_out.write("_atom_site_fract_z\n")

def read_if_positions(): # read initial or final positions 
    n_structures = 0
    for line in f_in:
        item = re.search(r"     number of atoms/cell      =",line)
        if item:
            n_atoms = int(re.search(r"\d+",line).group())
        
        item = re.search(r"lattice parameter \(alat\)\s+=\s+([\d\.]+)",line)
        if item:
            print("Found "+item.group())
            alat = float(re.search(r"-*\d+\.\d{4}",line).group())
            n_structures = n_structures+1
            
        item = re.search(r"crystal axes: \(cart\. coord\. in units of alat\)",line)
        if item:
            line = f_in.readline()
            a_vect = re.findall(r"-*\d+\.\d{6}",line)
            a_vect = [float(x)*alat*BOHR for x in a_vect]
            a_vect = np.asarray(a_vect)
            
            line = f_in.readline()
            b_vect = re.findall(r"-*\d+\.\d{6}",line)
            b_vect = [float(x)*alat*BOHR for x in b_vect]
            b_vect = np.asarray(b_vect)
            
            line = f_in.readline()
            c_vect = re.findall(r"-*\d+\.\d{6}",line)
            c_vect = [float(x)*alat*BOHR for x in c_vect]
            c_vect = np.asarray(c_vect)
                  
            a = np.linalg.norm(a_vect)
            b = np.linalg.norm(b_vect)
            c = np.linalg.norm(c_vect)
            alpha = math.degrees((cmath.acos(np.dot(b_vect, c_vect)/(b*c))).real)
            beta =  math.degrees((cmath.acos(np.dot(a_vect, c_vect)/(a*c))).real)
            gamma = math.degrees((cmath.acos(np.dot(a_vect, b_vect)/(a*b))).real)
            
            if n_structures == 1 : 
                data_title = "initial"
            else:
                data_title = "final"
            print_cif_header(a, b, c, alpha, beta, gamma, data_title)
            
        item = re.search(r"Cartesian axes",line)
        if item:
            line = f_in.readline() # skip line 
            line = f_in.readline() # skip line
                   
            for n in range(n_atoms):
                line = f_in.readline()
                el_name = re.search(r"[A-Z][a-z]{0,1}",line).group()
                atom = re.findall(r"-*\d+\.\d{7}",line)
                atom = [float(a)*alat*BOHR for a in atom]
                atom = np.asarray(atom)
                omega = np.dot(a_vect, np.cross(b_vect, c_vect))
                u = np.dot(atom, np.cross(b_vect,c_vect))/omega
                v = np.dot(atom, np.cross(c_vect,a_vect))/omega
                w = np.dot(atom, np.cross(a_vect,b_vect))/omega
                f_out.write("{} {:16.10f} {:16.10f} {:16.10f}\n".format(el_name,u,v,w))

def read_int_positions(): # read intermediate positions 
    n_structures = 0
    for line in f_in:
        item = re.search(r"     number of atoms/cell      =",line)
        if item:
            n_atoms = int(re.search(r"\d+",line).group())
        
        item = re.search(r"CELL_PARAMETERS \(alat=",line)
        if item:
            n_structures = n_structures+1
            print("Structure: "+str(n_structures))
            alat = float(re.search(r"\d+\.\d{8}",line).group())
            
            line = f_in.readline()
            a_vect = re.findall(r"-*\d+\.\d{9}",line)
            a_vect = [float(x)*alat*BOHR for x in a_vect]
            a_vect = np.asarray(a_vect)
            
            line = f_in.readline()
            b_vect = re.findall(r"-*\d+\.\d{9}",line)
            b_vect = [float(x)*alat*BOHR for x in b_vect]
            b_vect = np.asarray(b_vect)
            
            line = f_in.readline()
            c_vect = re.findall(r"-*\d+\.\d{9}",line)
            c_vect = [float(x)*alat*BOHR for x in c_vect]
            c_vect = np.asarray(c_vect)
                  
            a = np.linalg.norm(a_vect)
            b = np.linalg.norm(b_vect)
            c = np.linalg.norm(c_vect)
            alpha = math.degrees((cmath.acos(np.dot(b_vect, c_vect)/(b*c))).real)
            beta =  math.degrees((cmath.acos(np.dot(a_vect, c_vect)/(a*c))).real)
            gamma = math.degrees((cmath.acos(np.dot(a_vect, b_vect)/(a*b))).real)
            
            data_title = "structure_"+str(n_structures)
            print_cif_header(a, b, c, alpha, beta, gamma, data_title)
            
        item = re.search(r"ATOMIC_POSITIONS\s\(angstrom\)",line)
        if item:
            for n in range(n_atoms):
                line = f_in.readline()
                el_name = re.search(r"[A-Z][a-z]{0,1}",line).group()
                atom = re.findall(r"-*\d+\.\d{9}",line)
                atom = [float(a) for a in atom]
                atom = np.asarray(atom)
                omega = np.dot(a_vect, np.cross(b_vect, c_vect))
                u = np.dot(atom, np.cross(b_vect,c_vect))/omega
                v = np.dot(atom, np.cross(c_vect,a_vect))/omega
                w = np.dot(atom, np.cross(a_vect,b_vect))/omega
                f_out.write("{} {:16.10f} {:16.10f} {:16.10f}\n".format(el_name,u,v,w))

# constants 
BOHR = 0.52917721092 # Bohr radius, units of 10 nm 

cmd_args_parser = argparse.ArgumentParser()
cmd_args_parser.add_argument("-sf", help="starting and final structures",
                             action="store_true")
cmd_args_parser.add_argument("-int", help="intermediate structures",
                             action="store_true")
cmd_args_parser.add_argument("pw_out_file", help="pw.x output file")

cmd_args = cmd_args_parser.parse_args()

with open(cmd_args.pw_out_file, 'r') as f_in, \
open(cmd_args.pw_out_file+".cif",'w+') as f_out: 
    if cmd_args.sf: 
        read_if_positions()
    if cmd_args.int:
        read_int_positions()
         

    
