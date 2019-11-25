import bpy
import os
import sys
import json
import math

def main(elecFile=''):
	electrodes = open(elecFile)
	elecs = electrodes.readlines()
	f= open("guru99.txt","w+")
	f.write("This is line")


if __name__ == "__main__":
	main(sys.argv[1])
