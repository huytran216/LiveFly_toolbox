# **LiveFly** - A toolbox for the analysis of transcription dynamics in live Drosophila embryos

## Introduction

We present the LiveFly toolbox for the analysis of transcription dynamics in live Drosophila melanogaster embryos. The toolbox allows users to process two-color 3D confocal movies acquired using nuclei-labeling and the fluorescent RNA-tagging system described in the previous chapter and export the nuclei’s position as a function of time, their lineages and the intensity traces of the active loci.

The toolbox, which is tailored for the context of Drosophila early development, is semi-automatic, and requires minimal user intervention. It also includes a tool to combine data from multiple movies and visualize several features of the intensity traces and the expression pattern.

## System requirements

1. MATLAB 2013 or above, with the Image Processing Toolbox
2. Bioformat package for MATLAB: version 5.2.1 or higher.
3. The LiveFly toolbox: available for download from: https://github.com/huytran216/LiveFly_toolbox

## Input

1.	A time-lapse 3D movie captured with a confocal microscope. The compatible file formats are .czi or .lsm commonly used in time-lapse microscopy. The movie has two channels: one for the nuclei markers and one for MCP-GFP/PP7 proteins.
2.	Input parameters: time resolution and XYZ pixel resolution, embryo orientation and the coordinate of the imaging region relative to the embryo anterior/posterior poles.
3.	(Optional) List of nuclei IDs, time frames to be ignored for spot detection.

You can download a movie sample from: http://xfer.curie.fr/get/CDi5UvAjJZW/RAWMovie.tif (File size: 7.4 Gb)

The information on the movie can be found in: http://xfer.curie.fr/get/t5KwEYFFvrh/Movie%20information.txt (File size: 1 kb)

## Getting started

The toolbox has three modules:
1. The Nuclei segmentation module extracts from the 3D movie’s nuclei channel (marked by either fluorescent histones or NUP) the temporal position and lineage of each nucleus.
2. The Spot detection module extracts the intensity trace of active transcription loci in each nucleus.
3. The Visualizer module allows users to manage the data from multiple movies and to visualize the pattern of several features of the transcription dynamics.

Familiarize yourself with the tools by following the guide on Manual/Manual.docx

## License

Copyright (C) 2018 Huy Tran

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
