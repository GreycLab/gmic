---
title: "`G'MIC`: An Open-Source Self-Extending Framework for Image Processing"
tags:
  - image processing
  - image filtering
  - computer graphics
  - scripting language
  - user interfaces
  - creative coding
authors:
  - name: David Tschumperlé
    corresponding: true
    orcid: 0000-0003-3454-5079
    affiliation: 1
  - name: Sébastien Fourey
    affiliation: 1
    orcid: 0000-0001-9293-0771
  - name: Garry Osgood
    affiliation: 2
affiliations:
 - name: GREYC Lab (IMAGE Team), CNRS, Normandie Univ, UNICAEN, ENSICAEN, F-14000 Caen, France 
   index: 1
 - name: Independent contributor, New York City, US
   index: 2
date: 8 June 2023
bibliography: paper.bib
---

## Abstract

We present `G'MIC`, an open-source self-extending framework that defines an original, concise, scripting language for the writing of possibly complex
image processing operators and pipelines.
`G'MIC` also provides several user interfaces allowing for the manipulation of digital images, adapted to different levels of user expertise,
either from the command line, or as a _C/C`++`_ library, or as a user-friendly graphical plug-in that extends the capabilities of popular digital
image retouching applications, such as _GIMP_, _Krita_, _Photoshop_, _Affinity Photo_ and others.

## Keywords

Image Analysis, Processing and Filtering, Computer Graphics, Scripting Language, User Interfaces, Creative Coding.

# 1. Statement of need

## 1.1. Context

Tools not only shape work but set limits as to how far such work can go. Tools with limited capabilities or rigid behaviors not only circumscribe results but also constrain tool-users' awareness of what is possibile — they suppose from such inflexible tools only notions of what cannot be done.

Intrinsic to `G'MIC`'s design are means to map pipelines to commands, advancing the tool as a self-extending language and fortifying how users' conduct their work. Primal command pipelines may be further assembled into those having wider remits, these suitably named to bespeak their extended purposes and available for succeeding command prototyping.

`G'MIC` itself is based upon primitives of the _C`++`_ library [`CImg`](http://cimg.eu) [@cimg] that are broadly applicable to image processing work and which have been optimized for CPU performance. Most commands in the standard `G'MIC` distribution extend these primitives, using the aforementioned scheme. It is one which remains available to users, extending the language in line with their specific domains of expertise and making the language "their own".

`G'MIC` is distributed under the CeCILL free software licenses (GPL-compatible). The core language projects several user interfaces to convert, process or visualize generic *image datasets*. Allied with pipeline toolset, `G'MIC` embodies a highly flexible image model, ranging from 1D scalar signals to 3D+t sequences of multi-spectral volumetric images, hence including 2D color images.
This makes it a versatile tool for image processing, with a wide range of applications in research, industry and graphic design.

## 1.2. History and Motivation

The `G'MIC` project was started in mid-2008 by [David Tschumperlé](https://tschumperle.users.greyc.fr), a research scientist working in the IMAGE team of the _GREYC_, a public research laboratory affiliated with the CNRS institute in France.
David's area of research is the study and elaboration of image processing algorithms.

To that end, he first began developing [`CImg`](http://cimg.eu) [@cimg], beginning in 1999 and continuing to the present. `CImg` is  an open-source _C`++`_ library for generic image processing. Here, _generic_ implies a library that addresses structurally diverse imagery: photographs, multi-spectral images (e.g. from satellites), medical images (MRI, X-ray, tomography, etc.) and technical animations, among others.
The `CImg` library has therefore been designed to handle a wide variety of different image types, whether defined on 2D or 3D grids, or with any number of channels.
He still develops and maintains this library.

That said, `CImg` exhibits certain limitations for everyday research work:

1. When one simply wants to apply a predefined algorithm from `CImg` to an image, one needs to write a small, _C`++`_ program. Perhaps it is only a few lines long, but still it must be compiled and linked — and possibly debugged — before it can be executed. In the context of research work, such mechanics are just so many distractions. The idea of being able to run those algorithms directly from the command line is tempting.

2. Over time, a large number of these small, but purpose-specific, programs has accumulated. They solve specific problems but rarely see follow-on use. They are not broadly useful for integration into the `CImg` library and have become an unruly "collection" of specialized algorithms. By design, they cannot be easily distributed and are difficult to maintain (as opposed to a language having package managers, like Python).

These limitations motivated `G'MIC`'s development, beginning in 2008. Two design objectives came to the fore:

1. Enable _pipelines of image processing algorithms_ that may be directly invoked from the command line, without requiring compilation or linking steps.

2. Gather the implementation of specialized algorithms in a single location, facilitating their evolution, maintenance and distribution.

These objectives, in combination with a desire to write new image processing pipelines and algorithms in the most flexible and concise way possible, gave rise to the idea of _self-extension_.

It is well-known that research scientists are loathe to type; there are far more useful investments of time. In that light, there emerged a scheme to equate concise, shorthand words for pipelines — indeed, to _define_ new commands with pipelines. It is straightforward to see the advancement of this scheme, for such new commands can further define yet other commands. Such a mechanism also aligns well with a _write once, run everywhere_ doctrine. Improvements to a pipeline underlying a command propagates out to wherever that command is employed.

All these objectives led initially to the development of a specialized scripting language: the `G'MIC` language, and its associated interpreter, distributed as free software.

The first user interface created was `gmic`, the command line interface (_CLI_) tool that permits the execution of image processing code written in the `G'MIC` language directly from a shell. Other interfaces have followed since then, and will be detailed in Section 2.

## 1.4. Related Software

- **Command-line Interfaces:**

The command line interface `gmic` has been originally inspired by [_ImageMagick_](https://imagemagick.org/index.php) [@imagemagick] and [_GraphicsMagick_](http://www.graphicsmagick.org/) [@graphicsmagick], particularly the idea of being able to manipulate digital images from a shell. What all these projects have
in common is that they define distinct command languages, enabling the creation of image processing pipelines of varying complexity.

The main differences between `G'MIC` and _ImageMagick_/_GraphicsMagick_ are as follows:

1. The type of images processed is more diverse in `G'MIC`. Although _ImageMagick_ and _GraphicsMagick_ are capable to a certain extent of loading volumetric or hyperspectral images, the possibilities for processing these generic images is limited to the use of certain filters only
(on the other hand, _ImageMagick_ and _GraphicsMagick_ offer far more possibilities for converting image file formats, with format-specific encoding options).

2. The possibilities offered by the scripting languages associated with each project, for writing image processing pipelines, are more extensive in `G'MIC`. In particular, `G'MIC`'s scripting language makes it possible to write conditions, loops and multi-threaded pipelines, without having to resort to an external scripting language (such as `sh` or `bash`, which are typically used in conjunction with _ImageMagick_/_GraphicsMagick_). The richness of the `G'MIC` built-in scripting language (detailed in Section 3) ensures maximum portability of the developed pipelines between different architectures (_Linux_/_Windows_/_BSD_).

- **Image Filter Collections:**

There are also related software packages offering predefined filter sets to be applied to images. Popular examples are Mathmap [@mathmap], Filter Forge [@filterforge] and Pixelitor [@pixelitor].
While these software somehow allows the user to create its own pipeline of image processing filters, their use case is restricted to the provided graphical user interfaces, with quite limited scripting possibilities.

# 2. Framework Environment

## 2.1. Core Components

The current architecture of the `G'MIC` framework is depicted on Fig. 1. This corresponds to the current state of the framework (version **3.3.2**), at the time of writing.

![Overview of the `G'MIC` framework.](images/gmic_architecture.png)

The organization of this framework revolves around a central component: the **`G'MIC` scripting language interpreter** (in yellow). This interpreter uses the native functionalities of the **`CImg` library** (which is implemented in _C`++`_, in blue), but relies also on a set of commands, written in the `G'MIC` language themselves, constituting a **_standard library_ (`stdlib`)** for the framework (in green). The other components (in orange) stand for the different user interfaces provided by the framework.

More than 1000 distinct commands are currently implemented in the `stdlib`, covering a large portion of general image processing needs.
These commands are gathered by categories, and documented on the [reference pages](https://gmic.eu/reference/list_of_commands.html) of the project. The table below lists these categories, sorted by the respective number of commands they contain, and gives examples of typical commands found in each category:

| **Category**   | **# of commands** | **Examples of key commands** |
|:---|:---:|:---|
| Colors | 107 | `rgb2hsv`, `rgb2lab`, `retinex`, `sepia` |
| Filtering | 105 | `convolve`, `dilate`, `fft`, `sharpen` |
| Convenience Functions | 105 | `files`, `img2base64`, `strcapitalize` |
| 3D Meshes | 95 | `isosurface3d`, `rotate3d`, `torus3d` |
| Input / Output | 89 | `camera`, `echo`, `input`, `output`, `display` |
| Mathematical Operators | 58 | `add`, `argmax`, `cos`, `mul`, `sqrt` |
| Geometry Manipulation | 55 | `crop`, `resize`, `rotate`, `split` |
| Neural Networks | 56 | `nn_load`, `nn_conv2d`, `nn_maxpool2d` |
| Value Manipulation | 54 | `cut`, `equalize`, `normalize`, `map` |
| Interactive Commands | 47 | `demos`, `x_pacman`, `x_warp` |
| Features Extraction | 46 | `betti`, `histogram`, `label`, `skeleton` |
| Image Drawing | 41 | `ellipse`, `graph`, `line`, `polygon`, `text` |
| Artistic | 39 | `cartoon`, `cubism`, `polaroid`, `stencil` |
| Flow Control | 31 | `do`, `error`, `for`, `if`, `return`, `while` |
| Arrays, Tiles and Frames | 28 | `array`, `frame_xy`, `frame_blur` |
| Warpings | 24 | `deform`, `fisheye`, `twirl`, `warp` |
| Image Sequences and Videos | 20 | `animate`, `morph`, `apply_video` |
| Degradations | 13 | `cracks`, `pixelize`, `vignette` |
| Blending and Fading | 12 | `blend`, `fade_linear`, `fade_radial` |
| Matrix Computation | 11 | `dijkstra`, `eigen`, `invert`, `svd` |
| List Manipulation | 10 | `move`, `name`, `remove`, `reverse` |
| Other Commands | 3 | `debug`, `help`, `version` |

The `G'MIC` interpreter lets the user write and run custom programs using this predefined set of commands, for tasks as varied as writing new image filters, implementing generative algorithms or creating user interfaces for image manipulation.

## 2.2. User Interfaces

On top of the `G'MIC` interpreter are the user interfaces. Several types of user interface are implemented in the `G'MIC` framework, adapted to varying degrees of user's expertise.
Those interfaces are :

- **`gmic`**, a _command-line_ tool used to control the `G'MIC` interpreter from a terminal. It is actually one of the most powerful interface of the project, as it is able to manage all kind of image types (1D, 2D, 3D, multi-spectral, etc.). It can also open display windows for having basic user interaction when needed, typically for displaying images. Fig. 2 shows an example of use of `gmic` from a console, where a color image is imported, resized, blurred, converted to a 3D elevation mesh, and finally displayed in an interactive window.

![The command-line interface `gmic` in action.](images/gmic_cli2.png)

- **_G'MIC-Qt_** is a _Qt_-based [@qt] graphical interface intended to be used as a _plug-in_ for digital image retouching software, such as _GIMP_, _Krita_, _DigiKam_, _Photoshop_, _Affinity Photo_ and others, or as a _stand-alone_ program. This interface focuses on 2D color image processing. It proposes a set of filters (over 590 to date) to be applied to the user's input images. It features a fairly advanced system of dynamic user interface generation, based on the syntactic analysis of comment lines defined in `G'MIC` command files (thus including the image filters defined in the `stdlib`). All proposed filters are therefore written in the `G'MIC` language, with their parameter setting interface dynamically generated by the plug-in. It also features a filter update function, allowing to add/remove or correct existing filters without having to re-install new binaries of the software. Fig. 3 shows the _G'MIC-Qt_ interface applying an effect to a color image, here run from _GIMP_. Users of the plug-in are able to write their own `G'MIC` command files, in order to add new custom filters (with the corresponding _GUI_) into the plug-in.

![The _G'MIC-Qt_ plug-in in action.](images/gmic_qt_330.jpg)

- **`G'MIC` Online** is a website where a user can upload a color image and apply one of the _G'MIC-Qt_ filters on it (Fig.4).
It is a simple way to test the `G'MIC` filters and effects without having to install anything locally on the user's computer. It is written in CSS/Javascript and relies on the `gmic` _CLI_ tool on the server side to render the image filters.

![The _`G'MIC` Online_ website.](images/ui_gmicol_2.jpg)

- **`libgmic`** and **`libgmic`** are respectively _C`++`_ and _C_ libraries which allow the access to the `G'MIC` features directly from a _C/C`++`_ source code. They basically provide a simple _C/C`++`_ _API_ to run a `G'MIC` pipeline on a set of input images passed to the library.

- **_ZArt_** is a _Qt_-based graphical interface used mainly for demonstration purposes (Fig. 5), which applies `G'MIC` filters and effects on streamed webcam images in (almost) real-time.

![View of the _ZArt_ interface.](images/ui_zart.jpg)

- [**`gmic-py`**](https://pypi.org/project/gmic/) is a project for getting a _Python_ binding for `G'MIC` (still work-in-progress).
Its aim is to provide Python programmers with the full range of filters and image processing functions included in the `G'MIC` framework.

## 2.3. Visibility and Community

The `G'MIC` framework has been developed since 2008, mainly in the [IMAGE team](https://www.greyc.fr/equipes/image/) at the
[_GREYC_ laboratory](https://www.greyc.fr/), a French public research laboratory specialized in computer sciences.
The project web page is [https://gmic.eu](https://gmic.eu).
This website brings together a range of resources, from software download links to documentation and tutorial pages.

The core features of the `G'MIC` interpreter are developed by [David Tschumperlé](https://tschumperle.users.greyc.fr/),
the _G'MIC-Qt_ plug-in by [Sébastien Fourey](https://foureys.users.greyc.fr/), both being permanent researchers at _GREYC_.
The other contributors (for documentation, creation of new filters, or implementation of other user interfaces) can be found on
the [software's forum pages](https://discuss.pixls.us/c/software/gmic/10), hosted by [Pixls.Us](https://pixls.us/),
an association whose goal is to promote the use of open-source software dedicated to photography and image creation.
This forum is the place to go to get answers to questions about the software and chat with developers and users.

The `G'MIC` source code is available on these various github repositories:
[`gmic`](https://github.com/GreycLab/gmic/) (interpreter), [`gmic-qt`](https://github.com/c-koi/gmic-qt/) (plug-in) and
[`gmic-community`](https://github.com/GreycLab/gmic-community/) (external contributions, documentation).

Last but not least, the project provides regular updates on new developments on social networks such as
[Mastodon](https://piaille.fr/@gmic) and [Twitter](https://twitter.com/gmic_eu).

# 3. Examples of Research Work Conducted With `G'MIC`

To illustrate the high degree of genericity of the `G'MIC` framework, we list a selection of a few image processing research projects
carried out in the [IMAGE team](https://www.greyc.fr/equipes/image) at the [_GREYC_ laboratory](https://www.greyc.fr) (UMR CNRS 6072),
which have used `G'MIC` for algorithm development, prototyping, testing and result generation.

## 3.1. Patch-Based Image Inpainting

Between 2011 and 2015, the problem of Image _inpainting_ was studied by researchers Maxime Daisy, Pierre Buyssens, David Tschumperlé, Olivier Le Meur and Olivier Lézoray, and this led to the development of new, original image inpainting algorithms, described in detail
in a publication in _IEEE Transaction on Image Processing_ [@buyssens2015exemplar].
In _G'MIC-Qt_, the _Inpaint [Multi-Scale]_ and _Inpaint [Patch-Based]_ filters implement different patch-based image inpainting algorithms (Fig. 6).

![Patch-based image inpainting with `G'MIC`. Left: input image. Middle: user-defined mask. Right: inpainting result.](images/inpaint.png)

## 3.2. Color LUT Compression

3D _CLUTs_ (Color Look Up Tables) are popular digital models used in artistic image and video processing, for color grading,
simulation of analog films, and more generally for the description and application of generic non-parametric color transformations.
The relatively large size of these models leads to high data storage requirements when trying to distribute them on a large scale
(_e.g._ several hundred at the same time), typically as sets of `.png` HaldCLUTs or Adobe's `.cube` files.

In the context of `G'MIC`, it was important to be able to provide users with as many color transformation filters as possible,
without considerably increasing the size of the framework.
To that purpose, in 2018, researchers David Tschumperlé, Amal Mahboubi and Christine Porquet have proposed a dedicated method
for compressing _CLUTs_ [@tschumperle2020reconstruction], which delivers good performance (generally over 95% space saving) (Fig. 7).

![Principle of the `G'MIC` color _LUT_ compression algorithm. An input _CLUT_ (a) is analyzed and relevant color keypoints are deduced (b) and stored as a small image (c). A perceptual metric is used to ensure that the application of the compressed _CLUT_ on an image is visually similar to the application of the original one.](images/clut_compression2.png)

## 3.3. Semi-automatic Colorization of Line Arts

Colorizing line art drawings is a problem that illustrators are familiar with.
The question is how to colorize, with solid colors, an image originally made up of black or grayscale lines, on a white or transparent background.
The traditional tools available in image creation or retouching software (such as the well-known _Bucket Fill_) are not always well suited
because they do not take into account the specificities of the task, such as the fact that the lines of the drawing contain gaps which may be large,
or that the filling of colors should ideally be done under the lines (i.e. in a separate layer placed under the original grayscale anti-aliased drawing).

It was while discussing with David Revoy [@davidrevoy], an independent illustrator,
author of the webcomic _Pepper & Carrot_ [@pepperandcarrot],
that researchers Sébastien Fourey and David Tschumperlé came up with the idea of an algorithm which would make it possible to semi-automatically
generate a layer of colorization from an input line art.
The resulting algorithm analyzes the geometry of the line art contours and automatically deduces a reasonable flat-colored layer,
filled with colors which are randomly chosen, but that can be easily modified subsequently by the artist to give them the desired color
(for example, using the _Bucket Fill_ tool on the colorization layer thus generated).
This is the process illustrated in Fig. 9, where the input line art (_left_) is automatically colorized with random colors (_middle left_), which
are modified by the artist (_middle right_), before going to the illumination task (_right_).

![Principle of the semi-automatic line art colorization algorithm.](images/lineart_colorization.jpg)

The proposed method has been described in a conference paper,
published in the EUROGRAPHICS International Symposium on Vision, Modeling and Visualization, in 2018 [@fourey2018fast].

This filter offers three variants that can be used to colorize line arts. In addition to the mode of generating a layer containing random colors,
the filter offers a mode of intelligent extrapolation of color spots placed on a transparent layer above the original drawing,
as illustrated in Fig. 9.

![Color spot extrapolation for automatic lineart colorization.](images/lineart_cat.png)

Note that the colorization algorithm resulting from this research work was the subject of an external implementation in the GIMP software,
to enrich the "Bucket Fill" tool with a specialized "Line Art" mode for the colorization of line drawings [@jehan18].

## 3.4. Automatic Illumination of Flat-colored Drawings

In 2021, researchers David Tschumperlé, Amal Mahboubi and Christine Porquet have been interested in going one step further by
designing an original algorithm that tries to illuminate flat-colorized drawings, by automatically creating a light and shadow layer.
This method has been published and presented at the IEEE International Conference on Image Processing, in 2022 [@tschumperle2022automatic].

The main idea behind the algorithm is the analysis of the different silhouettes composing the drawing, such that plausible 3D elevation maps
are built. A Phong lightling model that relies on the corresponding normal maps is then applied to generate the illumination layer
(Fig. 10).

![Principle of the Shape Illumination Algorithm. Left: input image, middle-left: estimated 3D normal map. Right: two examples of different illuminations obtained with the Phong lighting model applied with different parameters.](images/illumination.png)

In the _G'MIC-Qt_ plug-in, this illumination algorithm can be applied via the filter **Illuminate 2D Shape**.

## 3.5. Patch-Based Image Style Transfer

Image stylization is a relatively recent processing application, having made its appearance in 2015 with the pioneering work of
[@gatys2015]. The problem consists in transforming an input image (usually a photograph) to give it a pictorial style close to
that of a second image (styled image), also specified by the user. This style image is, for example, the image of a famous painting,
a drawing or another photograph.

Classically, style transfer techniques are based on convolutional neural networks, but some lighter alternative methods exist.
Researchers Benjamin Samuth, David Tschumperlé and Julien Rabin have turned their attention to patch-based methods.
They were able to develop a multi-scale algorithm based solely on copying patches from the style image, to generate a coherent style transfer,
with low algorithmic cost.

This patch-based style transfer algorithm has been published and presented at the IEEE International Conference on Image Processing, in 2022 [@samuth2022patch].
It has been implemented in `G'MIC`, as command `stylize`, and its associated _G'MIC-Qt_ filter **Stylize** (Fig. 11).

![Examples of application of the `G'MIC` style transfer method. An input image (top left) is stylized according to a set of different style images (top row).](images/style_transfer.png)

## 3.6. Debanding of Astronomical Images

Regarding the use of G'MIC in research activities outside the _GREYC_ laboratory (in which it is developed) :
Let us mention that `G'MIC` is known to be used in the astronomy research community, in particular for processing images
from the [JWST (James Webb Space Telescope)](https://en.wikipedia.org/wiki/James_Webb_Space_Telescope),
which often exhibit frequency noise that materializes as transverse bands degrading the image quality.
One interesting algorithm for getting rid of band noise in `G'MIC` is implemented in the filter **Banding Denoise**.

For instance, the use of `G'MIC` is mentionned in the article [@ray2023outflows],
where images from protostar _HH211_ have been processed with it.
One of those made the cover of _Nature_ magazine (October 5, 2023, Volume 622 Issue 7981) (Fig. 12).

![Left: image of protostar _HH211_, partially processed with `G'MIC` (cover of Nature, courtesy of [Mark McCaughrean/ESA](https://mastodon.social/@markmccaughrean)). Right: an example of the effect of the `G'MIC` **Banding Denoise** filter on an image of the _IC4553_ galaxy (acquired by the JWST, courtesy of [Judy Schmidt](https://astrodon.social/@spacegeck)).](images/nature.png)

It is not possible to know the full range of uses for open-source software such as `G'MIC`, but it is reasonable to assume
it is used in other areas of research where image processing may be required.

# 5. Conclusions and Perspectives

In this article, we presented `G'MIC`, an open-source framework for digital image processing, developed since more than 15 years.
`G'MIC` defines its own scripting language to ease the design and application of image processing pipelines and
to allow the definition of new image processing filters and algorithms.
This mature software package can be used by a wide range of users : experts in image processing,
regular users looking for a tool to quickly retouch or generate images, or researchers working in fields where image processing can play a role to improve
acquired data.
The wide range of available user interfaces makes `G'MIC` a versatile tool for image processing.

With the advent of neural network-based image processing methods, `G'MIC`'s next major challenge will be to enable these new types
of processing to be used consistently within its scripting language.
Efforts are already underway to move in this direction.

# References
