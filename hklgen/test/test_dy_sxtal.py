import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
from copy import copy
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod
import sxtal_model as S
from pycrysfml import getSpaceGroup_crystalsys as xtalsys
import unittest
np.seterr(divide="ignore", invalid="ignore")

DATAPATH = os.path.dirname(os.path.abspath(__file__))
infoFile = os.path.join(DATAPATH,r"dy.cfl")
(spaceGroup, crystalCell, magAtomList, symmetry) = H.readMagInfo(infoFile)
spaceGroup, crystalCell, atomList = H.readInfo(infoFile)
wavelength = 1.703700
cell = crystalCell
refList = H.hklGen(spaceGroup, cell, np.float32(0.0), np.sin(179.5/2.)/wavelength, xtal=True)
magRefList = H.satelliteGen(cell, symmetry, np.sin(179.5/2.)/wavelength, hkls=refList)
def getMagStrFacts():
    return [float("%.3f" % I) for I in S.calcXtalIntensity(magRefList, magAtomList, symmetry, wavelength, cell, True)[0]]
def main():
    tt = [H.twoTheta(H.calcS(cell, ref.hkl),wavelength) for ref in magRefList]
    S.diffPatternXtal(infoFile=infoFile, cell=cell, scale=1, tt=tt, 
                      obsIntensity=np.zeros(len(magRefList)), wavelength=wavelength, basisSymmetry=symmetry, 
                      magAtomList=magAtomList, plot=True, residuals=False, nuclear=False, error=None, magnetic=True, 
                      info=True, base=0, refList=refList, extinctions=None)    
class TestMagStructureFactorsXtal(unittest.TestCase):
    def test_magstrFacts_dy_sxtal(self):
        """Test of Magnetic Structure Factor Calculations for single crystal Dy"""
        answers = [2.084, 2.039, 2.039, 1.995, 1.995, 0.829, 0.829, 0.829, 1.01, 1.01, 0.829, 1.01, 1.01, 1.168, 1.027, 1.027, 1.168, 1.027, 1.168, 1.027, 1.168, 1.871, 1.871, 1.144, 1.144, 1.259, 1.259, 1.259, 1.144, 1.144, 1.259, 1.231, 1.231, 1.231, 1.231, 1.795, 1.795, 1.278, 1.278, 1.278, 1.278, 1.306, 1.306, 1.306, 1.306, 0.557, 0.557, 0.557, 0.557, 1.293, 1.293, 1.365, 1.365, 1.293, 1.365, 1.365, 1.293, 0.661, 0.661, 0.661, 0.661, 1.373, 1.373, 1.317, 1.317, 1.317, 1.317, 1.373, 1.373, 0.739, 0.739, 0.739, 0.739, 1.33, 1.33, 1.33, 1.33, 1.619, 1.619, 1.32, 1.32, 1.32, 1.32, 0.882, 0.882, 0.882, 0.882, 1.526, 1.526, 0.789, 0.789, 0.646, 0.789, 0.789, 0.646, 0.646, 0.646, 0.929, 0.929, 0.929, 0.929, 1.323, 1.288, 1.323, 1.323, 1.288, 1.288, 1.323, 1.288, 0.816, 0.684, 0.684, 0.684, 0.684, 0.816, 0.816, 0.816, 0.971, 0.971, 0.971, 0.877, 0.971, 0.877, 0.877, 0.877, 0.714, 0.837, 0.837, 0.714, 0.837, 0.837, 0.714, 0.714, 0.892, 0.979, 0.892, 0.892, 0.979, 0.892, 0.979, 0.979, 1.256, 1.256, 1.256, 1.256, 1.275, 1.247, 1.247, 1.275, 1.247, 1.247, 1.275, 1.275, 0.983, 0.983, 0.902, 0.983, 0.983, 0.902, 0.902, 0.902, 0.775, 0.874, 0.775, 0.874, 0.775, 0.874, 0.775, 0.874, 0.967, 0.967, 0.967, 0.967, 1.208, 1.208, 1.208, 1.208, 0.502, 0.433, 0.502, 0.433, 0.502, 0.502, 0.433, 0.433, 1.333, 1.333, 0.981, 0.981, 0.981, 0.915, 0.915, 0.981, 0.915, 0.915, 0.884, 0.797, 0.884, 0.797, 0.797, 0.884, 0.884, 0.797, 0.529, 0.529, 0.464, 0.464, 0.529, 0.529, 0.464, 0.464, 0.552, 0.491, 0.491, 0.552, 0.491, 0.552, 0.552, 0.491, 0.961, 0.961, 0.961, 0.961, 0.973, 0.914, 0.973, 0.914, 0.914, 0.973, 0.973, 0.914, 1.133, 1.151, 1.133, 1.151, 1.151, 1.151, 1.133, 1.133, 1.239, 1.239, 0.6, 0.549, 0.6, 0.549, 0.549, 0.6, 0.6, 0.549, 0.828, 0.828, 0.828, 0.828, 0.813, 0.877, 0.813, 0.877, 0.813, 0.877, 0.813, 0.877, 0.827, 0.827, 0.827, 0.827, 0.621, 0.575, 0.621, 0.575, 0.621, 0.621, 0.575, 0.575, 1.09, 1.09, 1.09, 1.09, 0.825, 0.825, 0.825, 0.825, 0.891, 0.935, 0.891, 0.935, 0.891, 0.935, 0.891, 0.935, 1.082, 1.068, 1.068, 1.068, 1.082, 1.082, 1.082, 1.068, 0.806, 0.862, 0.806, 0.862, 0.806, 0.862, 0.806, 0.862, 0.911, 0.911, 0.911, 0.911, 0.813, 0.813, 0.813, 0.813, 0.869, 0.907, 0.907, 0.907, 0.869, 0.869, 0.869, 0.907, 1.026, 1.026, 1.026, 1.026, 0.608, 0.644, 0.644, 0.644, 0.608, 0.608, 0.608, 0.644, 0.802, 0.802, 0.802, 0.802, 0.874, 0.874, 0.874, 0.874, 1.055, 1.055, 0.55, 0.55, 0.45, 0.55, 0.45, 0.45, 0.45, 0.55, 0.335, 0.335, 0.335, 0.335, 0.571, 0.571, 0.571, 0.571, 0.658, 0.658, 0.658, 0.658, 0.554, 0.554, 0.554, 0.458, 0.458, 0.458, 0.458, 0.554, 0.348, 0.348, 0.348, 0.348, 0.614, 0.645, 0.645, 0.645, 0.614, 0.614, 0.645, 0.614, 0.574, 0.657, 0.574, 0.657, 0.574, 0.657, 0.574, 0.657, 0.466, 0.466, 0.558, 0.466, 0.466, 0.558, 0.558, 0.558, 0.36, 0.36, 0.36, 0.36, 0.768, 0.809, 0.809, 0.768, 0.768, 0.809, 0.809, 0.768, 0.655, 0.576, 0.655, 0.576, 0.576, 0.576, 0.655, 0.655, 0.94, 0.94, 0.94, 0.93, 0.93, 0.93, 0.93, 0.94, 0.768, 0.768, 0.768, 0.768, 0.835, 0.835, 0.807, 0.807, 0.835, 0.807, 0.835, 0.807, 0.97, 0.97, 0.482, 0.482, 0.563, 0.563, 0.563, 0.482, 0.563, 0.482, 0.387, 0.387, 0.387, 0.387, 0.647, 0.577, 0.647, 0.647, 0.577, 0.647, 0.577, 0.577, 0.774, 0.74, 0.774, 0.74, 0.774, 0.774, 0.74, 0.74, 0.893, 0.893, 0.893, 0.893, 0.746, 0.746, 0.746, 0.746, 0.563, 0.488, 0.563, 0.488, 0.488, 0.488, 0.563, 0.563, 0.319, 0.319, 0.385, 0.385, 0.385, 0.319, 0.385, 0.319, 0.401, 0.401, 0.401, 0.401, 0.64, 0.576, 0.64, 0.576, 0.576, 0.576, 0.64, 0.64, 0.771, 0.794, 0.794, 0.794, 0.771, 0.771, 0.771, 0.794, 0.328, 0.328, 0.328, 0.391, 0.328, 0.391, 0.391, 0.391, 0.604, 0.628, 0.604, 0.604, 0.628, 0.604, 0.628, 0.628, 0.87, 0.87, 0.863, 0.87, 0.863, 0.87, 0.863, 0.863, 0.783, 0.783, 0.783, 0.783, 0.579, 0.579, 0.579, 0.617, 0.617, 0.617, 0.579, 0.617, 0.335, 0.396, 0.335, 0.396, 0.396, 0.396, 0.335, 0.335, 0.576, 0.613, 0.613, 0.613, 0.576, 0.613, 0.576, 0.576, 0.574, 0.609, 0.574, 0.574, 0.574, 0.609, 0.609, 0.609, 0.828, 0.828, 0.828, 0.828, 0.555, 0.494, 0.494, 0.494, 0.494, 0.555, 0.555, 0.555, 0.408, 0.354, 0.408, 0.408, 0.354, 0.354, 0.408, 0.354, 0.42, 0.42, 0.42, 0.42, 0.591, 0.591, 0.611, 0.591, 0.611, 0.611, 0.611, 0.591, 0.618, 0.618, 0.564, 0.564, 0.618, 0.564, 0.618, 0.564, 0.735, 0.735, 0.735, 0.735, 0.69, 0.69, 0.69, 0.69, 0.564, 0.596, 0.596, 0.596, 0.564, 0.564, 0.564, 0.596, 0.413, 0.413, 0.362, 0.413, 0.362, 0.413, 0.362, 0.362, 0.491, 0.547, 0.491, 0.547, 0.547, 0.491, 0.547, 0.491, 0.425, 0.425, 0.425, 0.425, 0.671, 0.671, 0.671, 0.695, 0.695, 0.695, 0.671, 0.695, 0.602, 0.602, 0.602, 0.554, 0.602, 0.554, 0.554, 0.554, 0.811, 0.811, 0.557, 0.557, 0.586, 0.586, 0.586, 0.586, 0.557, 0.557, 0.706, 0.706, 0.706, 0.689, 0.689, 0.706, 0.689, 0.689, 0.659, 0.659, 0.659, 0.659, 0.375, 0.375, 0.417, 0.375, 0.417, 0.375, 0.417, 0.417, 0.737, 0.732, 0.737, 0.732, 0.732, 0.737, 0.737, 0.732, 0.633, 0.654, 0.654, 0.654, 0.633, 0.633, 0.654, 0.633, 0.536, 0.536, 0.561, 0.536, 0.561, 0.561, 0.561, 0.536, 0.549, 0.564, 0.564, 0.564, 0.564, 0.549, 0.549, 0.549, 0.74, 0.74, 0.521, 0.477, 0.521, 0.521, 0.477, 0.477, 0.477, 0.521, 0.647, 0.661, 0.661, 0.647, 0.661, 0.647, 0.661, 0.647, 0.378, 0.416, 0.378, 0.378, 0.378, 0.416, 0.416, 0.416, 0.255, 0.255, 0.255, 0.255, 0.231, 0.231, 0.231, 0.231, 0.424, 0.424, 0.424, 0.424, 0.704, 0.704, 0.704, 0.704, 0.563, 0.563, 0.526, 0.563, 0.563, 0.526, 0.526, 0.526, 0.236, 0.259, 0.259, 0.236, 0.259, 0.259, 0.236, 0.236, 0.522, 0.522, 0.545, 0.545, 0.545, 0.522, 0.522, 0.545, 0.263, 0.241, 0.241, 0.241, 0.263, 0.263, 0.241, 0.263, 0.636, 0.636, 0.636, 0.636, 0.476, 0.476, 0.476, 0.476, 0.672, 0.676, 0.672, 0.676, 0.672, 0.676, 0.672, 0.676, 0.465, 0.465, 0.465, 0.504, 0.504, 0.465, 0.504, 0.504, 0.418, 0.418, 0.418, 0.418, 0.473, 0.473, 0.473, 0.473, 0.537, 0.537, 0.524, 0.524, 0.537, 0.537, 0.524, 0.524, 0.541, 0.508, 0.541, 0.541, 0.541, 0.508, 0.508, 0.508, 0.591, 0.591, 0.591, 0.591, 0.355, 0.355, 0.417, 0.417, 0.417, 0.355, 0.417, 0.355, 0.252, 0.252, 0.272, 0.252, 0.272, 0.272, 0.252, 0.272, 0.356, 0.291, 0.291, 0.291, 0.356, 0.356, 0.291, 0.356, 0.469, 0.469, 0.469, 0.469, 0.355, 0.415, 0.415, 0.355, 0.415, 0.355, 0.415, 0.355, 0.355, 0.355, 0.292, 0.292, 0.355, 0.355, 0.292, 0.292, 0.646, 0.646, 0.646, 0.646, 0.355, 0.413, 0.355, 0.413, 0.355, 0.413, 0.355, 0.413, 0.258, 0.258, 0.258, 0.277, 0.277, 0.277, 0.258, 0.277, 0.375, 0.406, 0.406, 0.406, 0.375, 0.406, 0.375, 0.375, 0.294, 0.294, 0.294, 0.294, 0.355, 0.355, 0.355, 0.355, 0.458, 0.458, 0.458, 0.458, 0.589, 0.589, 0.589, 0.589, 0.489, 0.507, 0.507, 0.507, 0.507, 0.489, 0.489, 0.489, 0.555, 0.555, 0.555, 0.555, 0.569, 0.569, 0.569, 0.569, 0.556, 0.556, 0.556, 0.556, 0.353, 0.353, 0.406, 0.406, 0.406, 0.353, 0.406, 0.353, 0.353, 0.297, 0.297, 0.297, 0.353, 0.353, 0.297, 0.353, 0.451, 0.451, 0.451, 0.451, 0.398, 0.37, 0.398, 0.398, 0.398, 0.37, 0.37, 0.37, 0.562, 0.572, 0.572, 0.572, 0.562, 0.562, 0.562, 0.572, 0.351, 0.401, 0.401, 0.401, 0.351, 0.401, 0.351, 0.351, 0.283, 0.283, 0.283, 0.283, 0.267, 0.267, 0.267, 0.267, 0.37, 0.37, 0.37, 0.41, 0.37, 0.41, 0.41, 0.41, 0.351, 0.298, 0.351, 0.351, 0.351, 0.298, 0.298, 0.298, 0.611, 0.611, 0.463, 0.434, 0.434, 0.463, 0.434, 0.463, 0.463, 0.434, 0.487, 0.471, 0.487, 0.487, 0.471, 0.471, 0.471, 0.487, 0.407, 0.368, 0.407, 0.368, 0.368, 0.407, 0.368, 0.407, 0.398, 0.398, 0.398, 0.398, 0.491, 0.491, 0.466, 0.491, 0.466, 0.491, 0.466, 0.466, 0.264, 0.264, 0.216, 0.216, 0.264, 0.216, 0.264, 0.216, 0.366, 0.404, 0.366, 0.366, 0.404, 0.404, 0.366, 0.404, 0.265, 0.218, 0.265, 0.265, 0.218, 0.218, 0.218, 0.265, 0.285, 0.285, 0.285, 0.285, 0.27, 0.27, 0.27, 0.27, 0.516, 0.528, 0.528, 0.516, 0.516, 0.516, 0.528, 0.528, 0.432, 0.432, 0.432, 0.432, 0.469, 0.469, 0.478, 0.478, 0.478, 0.469, 0.478, 0.469, 0.265, 0.265, 0.265, 0.22, 0.22, 0.265, 0.22, 0.22, 0.561, 0.564, 0.564, 0.561, 0.564, 0.561, 0.561, 0.564, 0.344, 0.387, 0.387, 0.344, 0.344, 0.387, 0.344, 0.387, 0.53, 0.53, 0.521, 0.521, 0.521, 0.53, 0.53, 0.521, 0.36, 0.36, 0.395, 0.36, 0.395, 0.36, 0.395, 0.395, 0.343, 0.297, 0.297, 0.343, 0.297, 0.343, 0.343, 0.297, 0.416, 0.441, 0.416, 0.441, 0.416, 0.441, 0.416, 0.441, 0.42, 0.42, 0.42, 0.42, 0.384, 0.384, 0.384, 0.384, 0.443, 0.443, 0.443, 0.443, 0.465, 0.465, 0.465, 0.465, 0.225, 0.267, 0.225, 0.267, 0.225, 0.267, 0.225, 0.267, 0.356, 0.356, 0.356, 0.389, 0.389, 0.389, 0.389, 0.356, 0.54, 0.54, 0.54, 0.54, 0.555, 0.555, 0.339, 0.339, 0.339, 0.379, 0.339, 0.379, 0.379, 0.379, 0.352, 0.352, 0.374, 0.374, 0.374, 0.374, 0.352, 0.352, 0.338, 0.296, 0.338, 0.296, 0.296, 0.338, 0.338, 0.296, 0.267, 0.227, 0.267, 0.267, 0.227, 0.267, 0.227, 0.227, 0.486, 0.486, 0.486, 0.486, 0.448, 0.448, 0.44, 0.448, 0.44, 0.44, 0.44, 0.448, 0.498, 0.498, 0.498, 0.498, 0.429, 0.442, 0.429, 0.429, 0.442, 0.429, 0.442, 0.442, 0.27, 0.283, 0.27, 0.283, 0.27, 0.27, 0.283, 0.283, 0.17, 0.17, 0.17, 0.17, 0.514, 0.511, 0.511, 0.511, 0.514, 0.511, 0.514, 0.514, 0.172, 0.172, 0.172, 0.172, 0.344, 0.344, 0.373, 0.344, 0.344, 0.373, 0.373, 0.373, 0.174, 0.174, 0.174, 0.174, 0.36, 0.36, 0.36, 0.341, 0.36, 0.341, 0.341, 0.341, 0.393, 0.393, 0.393, 0.393, 0.231, 0.265, 0.265, 0.231, 0.265, 0.231, 0.231, 0.265, 0.268, 0.268, 0.28, 0.268, 0.28, 0.28, 0.28, 0.268, 0.493, 0.493, 0.493, 0.493, 0.453, 0.453, 0.453, 0.453, 0.358, 0.324, 0.358, 0.324, 0.358, 0.324, 0.358, 0.324, 0.407, 0.419, 0.407, 0.419, 0.407, 0.407, 0.419, 0.419, 0.337, 0.363, 0.363, 0.337, 0.337, 0.337, 0.363, 0.363, 0.324, 0.324, 0.324, 0.288, 0.288, 0.324, 0.288, 0.288, 0.343, 0.343, 0.343, 0.327, 0.343, 0.327, 0.327, 0.327, 0.179, 0.179, 0.179, 0.179, 0.442, 0.45, 0.45, 0.45, 0.442, 0.442, 0.442, 0.45, 0.34, 0.325, 0.34, 0.34, 0.325, 0.325, 0.34, 0.325, 0.457, 0.457, 0.457, 0.457, 0.263, 0.231, 0.231, 0.263, 0.263, 0.231, 0.231, 0.263, 0.394, 0.375, 0.394, 0.394, 0.375, 0.375, 0.375, 0.394, 0.378, 0.378, 0.378, 0.378, 0.351, 0.351, 0.351, 0.351, 0.337, 0.323, 0.323, 0.337, 0.337, 0.337, 0.323, 0.323, 0.395, 0.395, 0.395, 0.412, 0.412, 0.412, 0.395, 0.412, 0.181, 0.181, 0.181, 0.156, 0.181, 0.156, 0.156, 0.156, 0.181, 0.181, 0.181, 0.181, 0.45, 0.45, 0.444, 0.444, 0.45, 0.444, 0.444, 0.45, 0.182, 0.182, 0.158, 0.158, 0.158, 0.158, 0.182, 0.182, 0.346, 0.316, 0.346, 0.316, 0.346, 0.316, 0.346, 0.316, 0.315, 0.283, 0.283, 0.283, 0.315, 0.315, 0.315, 0.283, 0.159, 0.183, 0.159, 0.183, 0.159, 0.183, 0.183, 0.159, 0.329, 0.315, 0.315, 0.315, 0.315, 0.329, 0.329, 0.329, 0.341, 0.341, 0.341, 0.319, 0.319, 0.319, 0.319, 0.341, 0.311, 0.324, 0.324, 0.311, 0.324, 0.324, 0.311, 0.311, 0.388, 0.388, 0.383, 0.388, 0.383, 0.388, 0.383, 0.383, 0.163, 0.185, 0.185, 0.163, 0.185, 0.163, 0.163, 0.185, 0.454, 0.454, 0.371, 0.371, 0.371, 0.371, 0.354, 0.354, 0.354, 0.354, 0.185, 0.185, 0.185, 0.185, 0.259, 0.269, 0.259, 0.269, 0.269, 0.259, 0.259, 0.269, 0.415, 0.407, 0.415, 0.407, 0.407, 0.407, 0.415, 0.415, 0.333, 0.333, 0.333, 0.333, 0.371, 0.386, 0.371, 0.386, 0.371, 0.386, 0.386, 0.371, 0.229, 0.229, 0.229, 0.256, 0.229, 0.256, 0.256, 0.256, 0.313, 0.327, 0.327, 0.313, 0.313, 0.327, 0.313, 0.327, 0.165, 0.165, 0.186, 0.165, 0.186, 0.165, 0.186, 0.186, 0.414, 0.409, 0.409, 0.409, 0.414, 0.409, 0.414, 0.414, 0.186, 0.186, 0.186, 0.186, 0.363, 0.371, 0.371, 0.371, 0.363, 0.363, 0.363, 0.371, 0.309, 0.329, 0.329, 0.309, 0.309, 0.329, 0.309, 0.329, 0.346, 0.346, 0.346, 0.346, 0.424, 0.422, 0.422, 0.424, 0.424, 0.422, 0.424, 0.422, 0.215, 0.255, 0.255, 0.255, 0.215, 0.255, 0.215, 0.215, 0.31, 0.299, 0.31, 0.31, 0.299, 0.299, 0.31, 0.299, 0.226, 0.251, 0.226, 0.251, 0.251, 0.251, 0.226, 0.226, 0.389, 0.389, 0.389, 0.389, 0.319, 0.295, 0.295, 0.295, 0.319, 0.319, 0.319, 0.295, 0.215, 0.254, 0.215, 0.215, 0.254, 0.254, 0.215, 0.254, 0.261, 0.253, 0.253, 0.261, 0.261, 0.253, 0.261, 0.253, 0.294, 0.268, 0.268, 0.294, 0.268, 0.294, 0.294, 0.268, 0.221, 0.221, 0.18, 0.221, 0.18, 0.18, 0.18, 0.221, 0.214, 0.214, 0.253, 0.214, 0.253, 0.214, 0.253, 0.253, 0.355, 0.355, 0.36, 0.36, 0.36, 0.355, 0.36, 0.355, 0.408, 0.408, 0.408, 0.408, 0.22, 0.22, 0.22, 0.181, 0.181, 0.181, 0.22, 0.181, 0.168, 0.168, 0.168, 0.187, 0.187, 0.187, 0.187, 0.168, 0.298, 0.298, 0.298, 0.31, 0.31, 0.31, 0.298, 0.31, 0.291, 0.302, 0.302, 0.302, 0.291, 0.302, 0.291, 0.291, 0.219, 0.219, 0.219, 0.219, 0.181, 0.181, 0.181, 0.181, 0.411, 0.411, 0.329, 0.329, 0.329, 0.329, 0.382, 0.382, 0.382, 0.382, 0.341, 0.341, 0.348, 0.348, 0.341, 0.341, 0.348, 0.348, 0.248, 0.212, 0.248, 0.248, 0.248, 0.212, 0.212, 0.212, 0.169, 0.187, 0.169, 0.187, 0.169, 0.187, 0.169, 0.187, 0.283, 0.283, 0.305, 0.305, 0.283, 0.283, 0.305, 0.305, 0.186, 0.186, 0.186, 0.186, 0.259, 0.282, 0.259, 0.259, 0.282, 0.282, 0.259, 0.282, 0.258, 0.258, 0.258, 0.228, 0.228, 0.228, 0.258, 0.228, 0.217, 0.181, 0.181, 0.181, 0.181, 0.217, 0.217, 0.217, 0.211, 0.245, 0.245, 0.211, 0.211, 0.245, 0.211, 0.245, 0.256, 0.226, 0.256, 0.226, 0.226, 0.226, 0.256, 0.256, 0.385, 0.385, 0.384, 0.384, 0.385, 0.385, 0.384, 0.384, 0.359, 0.359, 0.359, 0.359, 0.286, 0.302, 0.286, 0.286, 0.302, 0.286, 0.302, 0.302, 0.225, 0.225, 0.254, 0.225, 0.225, 0.254, 0.254, 0.254, 0.215, 0.215, 0.18, 0.18, 0.215, 0.18, 0.215, 0.18, 0.324, 0.312, 0.312, 0.312, 0.324, 0.324, 0.324, 0.312, 0.295, 0.295, 0.295, 0.295, 0.324, 0.335, 0.335, 0.324, 0.324, 0.324, 0.335, 0.335, 0.218, 0.238, 0.218, 0.218, 0.218, 0.238, 0.238, 0.238, 0.185, 0.185, 0.185, 0.185, 0.371, 0.371, 0.371, 0.371, 0.284, 0.284, 0.275, 0.275, 0.284, 0.275, 0.275, 0.284, 0.17, 0.17, 0.17, 0.17, 0.139, 0.139, 0.139, 0.139, 0.222, 0.222, 0.222, 0.249, 0.222, 0.249, 0.249, 0.249, 0.348, 0.348, 0.348, 0.343, 0.343, 0.348, 0.343, 0.343, 0.254, 0.254, 0.254, 0.254, 0.17, 0.139, 0.139, 0.17, 0.17, 0.17, 0.139, 0.139, 0.237, 0.206, 0.237, 0.237, 0.237, 0.206, 0.206, 0.206, 0.243, 0.243, 0.243, 0.243, 0.236, 0.236, 0.236, 0.236, 0.252, 0.252, 0.252, 0.252, 0.349, 0.349, 0.349, 0.349, 0.169, 0.169, 0.184, 0.184, 0.184, 0.169, 0.184, 0.169, 0.17, 0.17, 0.14, 0.14, 0.14, 0.14, 0.17, 0.17, 0.347, 0.347, 0.347, 0.347, 0.344, 0.344, 0.344, 0.344, 0.289, 0.289, 0.274, 0.289, 0.274, 0.274, 0.289, 0.274, 0.245, 0.245, 0.219, 0.245, 0.219, 0.219, 0.245, 0.219, 0.25, 0.25, 0.25, 0.25, 0.179, 0.21, 0.179, 0.21, 0.179, 0.21, 0.179, 0.21, 0.266, 0.274, 0.274, 0.274, 0.266, 0.274, 0.266, 0.266, 0.212, 0.231, 0.231, 0.231, 0.212, 0.212, 0.212, 0.231, 0.122, 0.122, 0.113, 0.113, 0.122, 0.113, 0.113, 0.122, 0.233, 0.233, 0.204, 0.204, 0.233, 0.204, 0.204, 0.233, 0.295, 0.295, 0.295, 0.295, 0.114, 0.122, 0.114, 0.122, 0.122, 0.114, 0.114, 0.122, 0.291, 0.291, 0.302, 0.302, 0.302, 0.302, 0.291, 0.291, 0.169, 0.169, 0.141, 0.141, 0.141, 0.169, 0.141, 0.169, 0.276, 0.266, 0.266, 0.276, 0.276, 0.276, 0.266, 0.266, 0.277, 0.277, 0.277, 0.277, 0.302, 0.302, 0.311, 0.311, 0.311, 0.302, 0.311, 0.302, 0.245, 0.245, 0.245, 0.245, 0.177, 0.207, 0.207, 0.177, 0.177, 0.177, 0.207, 0.207, 0.123, 0.115, 0.115, 0.115, 0.123, 0.123, 0.115, 0.123, 0.182, 0.168, 0.168, 0.182, 0.182, 0.168, 0.182, 0.168, 0.258, 0.275, 0.258, 0.258, 0.258, 0.275, 0.275, 0.275, 0.257, 0.257, 0.239, 0.239, 0.239, 0.257, 0.257, 0.239, 0.303, 0.303, 0.307, 0.303, 0.307, 0.307, 0.307, 0.303, 0.141, 0.169, 0.141, 0.141, 0.141, 0.169, 0.169, 0.169, 0.304, 0.304, 0.298, 0.304, 0.298, 0.298, 0.304, 0.298, 0.233, 0.233, 0.227, 0.227, 0.227, 0.233, 0.233, 0.227, 0.236, 0.236, 0.236, 0.213, 0.213, 0.213, 0.213, 0.236, 0.241, 0.241, 0.241, 0.241, 0.18, 0.18, 0.18, 0.18, 0.319, 0.319, 0.315, 0.315, 0.315, 0.315, 0.319, 0.319, 0.124, 0.124, 0.124, 0.117, 0.117, 0.124, 0.117, 0.117, 0.211, 0.228, 0.228, 0.211, 0.228, 0.211, 0.211, 0.228, 0.336, 0.336, 0.318, 0.318, 0.315, 0.315, 0.315, 0.318, 0.318, 0.315, 0.231, 0.231, 0.209, 0.209, 0.231, 0.209, 0.231, 0.209, 0.278, 0.278, 0.278, 0.278, 0.227, 0.227, 0.227, 0.21, 0.21, 0.21, 0.227, 0.21, 0.125, 0.125, 0.118, 0.118, 0.118, 0.125, 0.125, 0.118, 0.221, 0.196, 0.221, 0.221, 0.196, 0.196, 0.196, 0.221, 0.208, 0.225, 0.208, 0.208, 0.225, 0.225, 0.208, 0.225, 0.304, 0.304, 0.304, 0.304, 0.261, 0.245, 0.245, 0.261, 0.245, 0.245, 0.261, 0.261, 0.166, 0.141, 0.141, 0.141, 0.166, 0.166, 0.166, 0.141, 0.25, 0.25, 0.25, 0.259, 0.259, 0.259, 0.259, 0.25, 0.228, 0.228, 0.228, 0.244, 0.244, 0.244, 0.244, 0.228, 0.253, 0.246, 0.246, 0.253, 0.246, 0.253, 0.246, 0.253, 0.231, 0.231, 0.231, 0.231, 0.173, 0.173, 0.198, 0.198, 0.198, 0.173, 0.198, 0.173, 0.177, 0.177, 0.177, 0.177, 0.248, 0.248, 0.26, 0.26, 0.26, 0.248, 0.26, 0.248, 0.316, 0.315, 0.315, 0.315, 0.316, 0.315, 0.316, 0.316, 0.283, 0.278, 0.283, 0.278, 0.278, 0.278, 0.283, 0.283, 0.283, 0.283, 0.28, 0.283, 0.28, 0.283, 0.28, 0.28, 0.22, 0.204, 0.22, 0.22, 0.204, 0.22, 0.204, 0.204, 0.199, 0.199, 0.199, 0.213, 0.199, 0.213, 0.213, 0.213, 0.164, 0.175, 0.164, 0.175, 0.164, 0.175, 0.175, 0.164, 0.164, 0.164, 0.141, 0.141, 0.141, 0.141, 0.164, 0.164, 0.192, 0.192, 0.192, 0.215, 0.215, 0.215, 0.192, 0.215, 0.226, 0.226, 0.226, 0.226, 0.119, 0.126, 0.126, 0.126, 0.119, 0.119, 0.119, 0.126, 0.202, 0.202, 0.202, 0.217, 0.217, 0.217, 0.202, 0.217, 0.305, 0.305, 0.305, 0.305, 0.194, 0.17, 0.17, 0.194, 0.194, 0.17, 0.194, 0.17, 0.122, 0.122, 0.103, 0.103, 0.103, 0.122, 0.103, 0.122, 0.219, 0.219, 0.199, 0.199, 0.219, 0.199, 0.199, 0.219, 0.242, 0.235, 0.242, 0.235, 0.242, 0.235, 0.235, 0.242, 0.103, 0.122, 0.103, 0.103, 0.122, 0.103, 0.122, 0.122, 0.126, 0.126, 0.126, 0.12, 0.12, 0.12, 0.126, 0.12, 0.303, 0.303, 0.289, 0.289, 0.289, 0.289, 0.122, 0.104, 0.104, 0.122, 0.122, 0.122, 0.104, 0.104, 0.246, 0.235, 0.246, 0.235, 0.235, 0.246, 0.235, 0.246, 0.28, 0.28, 0.28, 0.28, 0.161, 0.171, 0.171, 0.161, 0.171, 0.171, 0.161, 0.161, 0.252, 0.252, 0.26, 0.252, 0.26, 0.26, 0.26, 0.252, 0.206, 0.206, 0.211, 0.206, 0.211, 0.211, 0.211, 0.206, 0.241, 0.241, 0.241, 0.241, 0.191, 0.204, 0.204, 0.191, 0.191, 0.204, 0.204, 0.191, 0.266, 0.259, 0.259, 0.266, 0.266, 0.266, 0.259, 0.259, 0.194, 0.194, 0.194, 0.194, 0.212, 0.212, 0.212, 0.212, 0.195, 0.195, 0.195, 0.208, 0.208, 0.208, 0.208, 0.195, 0.139, 0.139, 0.16, 0.139, 0.16, 0.16, 0.16, 0.139, 0.122, 0.122, 0.105, 0.105, 0.105, 0.122, 0.122, 0.105, 0.213, 0.213, 0.213, 0.213, 0.286, 0.286, 0.286, 0.286, 0.286, 0.286, 0.286, 0.286, 0.122, 0.122, 0.122, 0.122, 0.105, 0.105, 0.105, 0.105, 0.191, 0.204, 0.191, 0.204, 0.191, 0.191, 0.204, 0.204]
        self.assertEqual(answers, getMagStrFacts())
if __name__ == "__main__":
    # program run normally
    main()