<a href="https://gmic.eu">![Logo](https://gmic.eu/img/logo4.jpg)</a>

#### A Full-Featured Open-Source Framework for Image Processing
##### https://gmic.eu

<a href="https://gmic.eu">![Logo](https://gmic.eu/img/logos2.jpg)</a>

-----------------------

## Overview

**G'MIC** is a full-featured open-source framework for **[digital image processing](https://en.wikipedia.org/wiki/Image_processing)**, distributed under the
**[CeCILL](http://cecill.info/index.en.html)** free software licenses (LGPL-like and/or GPL-compatible).
It provides several different **[user interfaces](https://en.wikipedia.org/wiki/User_interface)** for converting, manipulating, filtering and visualizing **generic image datasets**,
ranging from _1D scalar signals_ to _3D+t sequences of multi-spectral volumetric images_, which naturally includes _2D color images_.
These user interfaces are:

- **1** - **[gmic](https://gmic.eu/reference/)**, a command-line interface for using G'MIC's image processing features from a shell.
In this context, G'MIC can be viewed as a friendly companion to the [ImageMagick](http://www.imagemagick.org) or [GraphicsMagick](http://www.graphicsmagick.org) software suites.

<a href="https://gmic.eu/img/gmic_cli.jpg">![cli_interface](https://gmic.eu/img/gmic_cli_thumb.jpg)</a>

- **2** - **[libgmic](https://gmic.eu/libgmic.html)**, a small, portable, thread-safe, multi-threaded _C++_ image processing library designed to be linked with third-party applications.
[Its simple API](https://gmic.eu/libgmic.html) allows developers to integrate all **G'MIC** features into their own software without much effort (a _C_ API is available as well).

<a href="https://gmic.eu/img/gmic_libgmic.jpg">![libgmic_interface](https://gmic.eu/img/gmic_libgmic_thumb.jpg)</a>

- **3** - **[G'MIC-Qt](https://github.com/GreycLab/gmic-qt)**, a plug-in designed to bring **G'MIC**'s capabilities to image editing applications like
  [GIMP](http://www.gimp.org),
  [Krita](https://www.krita.org),
  [Photoshop](https://en.wikipedia.org/wiki/Adobe_Photoshop),
  [Affinity Photo](https://en.wikipedia.org/wiki/Affinity_Photo),
  [PaintShop Pro](https://en.wikipedia.org/wiki/PaintShop_Pro),
  [PhotoLine](https://en.wikipedia.org/wiki/PhotoLine) and
  [Paint.NET](https://www.getpaint.net/).
  More than **600 filters** are already available, sorted by category (_Artistic, Black & white, Colors, Contours, Deformations, Degradations, Details, Film emulation, Frames, Layers,
Light & shadows, Patterns, Rendering, Repair, Sequences_, etc.).

<a href="https://gmic.eu/img/gmic_gimp.jpg">![gimp_interface](https://gmic.eu/img/gmic_gimp_thumb.jpg)</a>

- **4** - **[G'MIC Online](https://gmicol.greyc.fr)**, a web service that lets users apply image processing algorithms to their images, directly from a web browser.

<a href="https://gmic.eu/img/gmic_gmicol.jpg">![gmicol_interface](https://gmic.eu/img/gmic_gmicol_thumb.jpg)</a>

- **5** - **[ZArt](https://www.youtube.com/watch?v=k1l3RdvwHeM)**, a Qt-based interface for real-time processing of video streams from webcams or video files.

<a href="https://gmic.eu/img/gmic_zart.jpg">![zart_interface](https://gmic.eu/img/gmic_zart_thumb.jpg)</a>

-----------------------

Other open-source projects are known to use some of the **G'MIC** features:

- **[EKD](http://ekd.tuxfamily.org)**, free software dedicated to post-production processing for videos and images.
- **[Flowblade](https://github.com/jliljebl/flowblade)**, a multitrack non-linear video editor for Linux released under GPL 3 license.
- **[Photoflow](http://aferrero2707.github.io/PhotoFlow/)**, a fully non-destructive photo retouching program.
- **[Veejay](http://veejayhq.net/)**, a visual instrument and realtime video sampler.

-----------------------

## Packaging Status

[![Packaging status](https://repology.org/badge/tiny-repos/gmic.svg)](https://repology.org/project/gmic/versions)
[![latest packaged version(s)](https://repology.org/badge/latest-versions/gmic.svg)](https://repology.org/project/gmic/versions)
<details>
  <summary>Expand to see the packaging status of G'MIC...</summary>

[![Packaging status](https://repology.org/badge/vertical-allrepos/gmic.svg?header=gmic)](https://repology.org/project/gmic/versions)

</details>

-----------------------

## Testing the software

You can verify that G'MIC is working correctly by running these two commands,
which execute all the commands and filters available in the framework:
~~~
$ gmic input_text https://gmic.eu/gmic_stdlib.\$_version parse_cli images
~~~
and
~~~
$ gmic input_text https://gmic.eu/gmic_stdlib.\$_version parse_gui images
~~~
Please make sure you run these commands in an empty folder, as they generate a lot of new image files.

-----------------------

## Referencing the software

To cite G'MIC in your own publications, please cite the JOSS article, published in January 2025:

- [D. Tschumperlé](https://tschumperle.users.greyc.fr),
  [S. Fourey](https://foureys.users.greyc.fr) and G. Osgood (2025).
  [G'MIC: An Open-Source Self-Extending Framework for Image Processing](https://doi.org/10.21105/joss.06618),
  Journal of Open Source Software, 10(105), 6618, https://doi.org/10.21105/joss.06618

- Corresponding _BibTeX_ entry:
~~~
@article{Tschumperlé2025,
  doi = {10.21105/joss.06618},
  url = {https://doi.org/10.21105/joss.06618},
  year = {2025},
  publisher = {The Open Journal},
  volume = {10},
  number = {105},
  pages = {6618},
  author = {David Tschumperlé and Sébastien Fourey and Garry Osgood},
  title = {G'MIC: An Open-Source Self-Extending Framework for Image Processing},
  journal = {Journal of Open Source Software}
}
~~~

-----------------------
