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

# 1. Statement of Need

## 1.1. Context

Intrinsic to `G'MIC`'s design are means to map image processing pipelines to commands, advancing the tool as a self-extending language
and fortifying how users' conduct their work. Primal command pipelines may be further assembled into those having wider remits, these suitably
named to bespeak their extended purposes and available for succeeding command prototyping.

`G'MIC` is distributed under the CeCILL free software licenses (GPL-compatible). The core language projects several user interfaces to convert,
process or visualize generic *image datasets*. Allied with pipeline toolset, `G'MIC` embodies a highly flexible image model,
ranging from 1D scalar signals to 3D+t sequences of multi-spectral volumetric images, hence including 2D color images.
This makes it a versatile tool for image processing, with a wide range of applications in research, industry and graphic design.

## 1.2. History and Motivation

The `G'MIC` project was started in mid-2008 by [David Tschumperlé](https://tschumperle.users.greyc.fr), a research scientist working in the IMAGE team of the _GREYC_, a public research laboratory affiliated with the CNRS institute in France.
David's area of research is the study and elaboration of image processing algorithms.

To that end, he first began developing [`CImg`](http://cimg.eu) [@cimg], beginning in 1999 and continuing to the present. `CImg` is  an open-source _C`++`_ library for generic image processing. Here, _generic_ implies a library that addresses structurally diverse imagery: photographs, multi-spectral images (e.g. from satellites), medical images (MRI, X-ray, tomography, etc.) and technical animations, among others.

That said, `CImg` exhibits certain limitations for everyday research work:

1. When one simply wants to apply a predefined algorithm from `CImg` to an image, one needs to write a small, _C`++`_ program. Perhaps it is only a few lines long, but still it must be compiled and linked — and possibly debugged — before it can be executed. In the context of research work, such mechanics are just so many distractions. The idea of being able to run those algorithms directly from the command line is tempting.

2. Over time, a large number of these small, but purpose-specific, programs has accumulated. They solve specific problems but rarely see follow-on use. They are not broadly useful for integration into the `CImg` library and have become an unruly "collection" of specialized algorithms. By design, they cannot be easily distributed and are difficult to maintain (as opposed to a language having package managers, like Python).

These limitations motivated `G'MIC`'s development, beginning in 2008. Two design objectives came to the fore:

1. Enable _pipelines of image processing algorithms_ that may be directly invoked from the command line, without requiring compilation or linking steps.

2. Gather the implementation of specialized algorithms in a single location, facilitating their evolution, maintenance and distribution.

These objectives, in combination with a desire to write new image processing pipelines and algorithms in the most flexible and concise way possible, gave rise to the idea of _self-extension_. All these objectives led initially to the development of a specialized scripting language: the `G'MIC` language, and its associated interpreter, distributed as free software.

## 1.4. Related Software

- **Command-line Interfaces:**

The command line interface `gmic` has been originally inspired by [_ImageMagick_](https://imagemagick.org/index.php) [@imagemagick] and [_GraphicsMagick_](http://www.graphicsmagick.org/) [@graphicsmagick], particularly the idea of being able to manipulate digital images from a shell. What all these projects have
in common is that they define distinct command languages, enabling the creation of image processing pipelines of varying complexity.

The main differences between `G'MIC` and _ImageMagick_/_GraphicsMagick_ are as follows:

1. The type of images processed is more diverse in `G'MIC`. Although _ImageMagick_ and _GraphicsMagick_ are capable to a certain extent of loading volumetric or hyperspectral images, the possibilities for processing these generic images is limited to the use of certain filters only

2. The possibilities offered by the scripting languages associated with each project, for writing image processing pipelines, are more extensive in `G'MIC`. In particular, `G'MIC`'s scripting language makes it possible to write conditions, loops and multi-threaded pipelines, without having to resort to an external scripting language (such as `sh` or `bash`, which are typically used in conjunction with _ImageMagick_/_GraphicsMagick_).

- **Image Filter Collections:**

There are also related software packages offering predefined filter sets to be applied to images. Popular examples are Mathmap [@mathmap], Filter Forge [@filterforge] and Pixelitor [@pixelitor].
While these software somehow allows the user to create its own pipeline of image processing filters, their use case is restricted to the provided graphical user interfaces, with quite limited scripting possibilities.

# 2. Framework Environment

## 2.1. Core Components

The current architecture of the `G'MIC` framework is depicted on Fig. 1. This corresponds to the current state of the framework (version **3.3.2**), at the time of writing.

![Overview of the `G'MIC` framework.](images/gmic_architecture.png)

The organization of this framework revolves around a central component: the **`G'MIC` scripting language interpreter** (in yellow). This interpreter uses the native functionalities of the **`CImg` library** (which is implemented in _C`++`_, in blue), but relies also on a set of commands, written in the `G'MIC` language themselves, constituting a **_standard library_ (`stdlib`)** for the framework (in green). The other components (in orange) stand for the different user interfaces provided by the framework. More than 1000 distinct commands are currently implemented in the `stdlib`, covering a large portion of general image processing needs.

The `G'MIC` interpreter lets the user write and run custom programs using this predefined set of commands, for tasks as varied as writing new image filters, implementing generative algorithms or creating user interfaces for image manipulation.

## 2.2. User Interfaces

On top of the `G'MIC` interpreter are the user interfaces. Several types of UI are implemented, adapted to varying degrees of user's expertise:

- **`gmic`**, a _command-line_ tool used to control the `G'MIC` interpreter from a terminal (Fig. 2).

![The command-line interface `gmic` in action.](images/gmic_cli2.png)

- **_G'MIC-Qt_** is a _Qt_-based [@qt] graphical interface intended to be used as a _plug-in_ for digital image retouching software, such as _GIMP_, _Krita_, _DigiKam_, _Photoshop_, _Affinity Photo_ and others, or as a _stand-alone_ program (Fig. 3).

![The _G'MIC-Qt_ plug-in in action.](images/gmic_qt_330.jpg)

- **`G'MIC` Online** is a website where a user can upload a color image and apply one of the _G'MIC-Qt_ filters on it.

- **`libgmic`** and **`libgmic`** are respectively _C`++`_ and _C_ libraries which basically provide a simple _C/C`++`_ _API_ to run a `G'MIC` pipeline on a set of input images.

- **_ZArt_** is a _Qt_-based graphical interface used mainly for demonstration purposes, which applies `G'MIC` filters and effects on streamed webcam images in (almost) real-time.

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
carried out, which have used `G'MIC` for algorithm development, prototyping, testing and result generation.

## 3.1. Patch-Based Image Inpainting

Between 2011 and 2015, the problem of Image _inpainting_ was studied by researchers Maxime Daisy, Pierre Buyssens, David Tschumperlé, Olivier Le Meur and Olivier Lézoray, in [@buyssens2015exemplar]. In `G'MIC`, several algorithms are implemented for image inpainting (Fig. 4).

![Patch-based image inpainting with `G'MIC`. Left: input image. Middle: user-defined mask. Right: inpainting result.](images/inpaint.png)

## 3.2. Color LUT Compression

3D _CLUTs_ (Color Look Up Tables) are popular digital models used in artistic image and video processing, for the description
of generic non-parametric color transformations.
The relatively large size of these models leads to high data storage requirements when trying to distribute them on a large scale
(_e.g._ several hundred at the same time), typically as sets of `.png` HaldCLUTs or Adobe's `.cube` files.
For storage purposes, researchers David Tschumperlé, Amal Mahboubi and Christine Porquet have proposed an original algorithm
for _CLUT_ compression [@tschumperle2020reconstruction], which allows `G'MIC` to provide more than 1100 _CLUTs_ in approximatively only 4MiB of storage (Fig. 5).

![Principle of the `G'MIC` color _LUT_ compression algorithm. An input _CLUT_ (a) is analyzed and relevant color keypoints are deduced (b) and stored as a small image (c). A perceptual metric is used to ensure that the application of the compressed _CLUT_ on an image is visually similar to the application of the original one.](images/clut_compression2.png)

## 3.3. Semi-automatic Colorization of Line Arts

Colorizing line art drawings is a problem that illustrators are familiar with.
Traditional tools available in image creation or retouching software (such as the well-known _Bucket Fill_) are not always well suited
because they do not take into account the anti-aliased nature of the lines, or the gaps that may be present in line drawings.
It was while discussing with David Revoy [@davidrevoy], an independent illustrator, author of the webcomic _Pepper & Carrot_ [@pepperandcarrot],
that researchers Sébastien Fourey and David Tschumperlé came up with the idea of an algorithm which would make it possible to semi-automatically
generate a layer of colorization from an input line art [@fourey2018fast].
The resulting _"Smart Coloring"_ algorithm, now implemented in `G'MIC`, analyzes the geometry of the line art contours and automatically deduces a reasonable flat-colored layer, from a user-defined layer that only contains a few color strokes (Fig. 6).

![Color spot extrapolation for automatic lineart colorization.](images/lineart_cat.png)

Note that this colorization algorithm has been also implemented natively in the GIMP software, to enrich the "Bucket Fill" tool with a specialized "Line Art" mode for the colorization of line drawings [@jehan18].

## 3.4. Automatic Illumination of Flat-colored Drawings

In a similar vein, researchers David Tschumperlé, Amal Mahboubi and Christine Porquet have been interested in going one step further by
designing an original algorithm that illuminates flat-colorized drawings, by automatically creating a light and shadow layer
[@tschumperle2022automatic] (Fig. 7).

![Principle of the Shape Illumination Algorithm. Left: input image, middle-left: estimated 3D normal map. Right: two examples of different illuminations obtained with the Phong lighting model applied with different parameters.](images/illumination.png)

## 3.5. Patch-Based Image Style Transfer

Image stylization consists in transforming an input image to give it a pictorial style close to that of a second image (style image).
In 2022, researchers Benjamin Samuth, David Tschumperlé and Julien Rabin have turned their attention to patch-based methods and
were able to develop a patch-based multi-scale algorithm, with a low algorithmic cost [@samuth2022patch].
It has been implemented in the `G'MIC` framework, as a new command `stylize` (Fig. 8).

![Examples of application of the `G'MIC` style transfer method. An input image (top left) is stylized according to a set of different style images (top row).](images/style_transfer2.png)

## 3.6. Debanding of Astronomical Images

Let us finally mention that `G'MIC` is known to be used in the astronomy research community, in particular for processing images
from the [JWST (James Webb Space Telescope)](https://en.wikipedia.org/wiki/James_Webb_Space_Telescope),
which often exhibit frequency noise that materializes as transverse bands degrading the image quality
(quite well removed with `G'MIC` filter **Banding Denoise**).
`G'MIC` has been mentionned in the article [@ray2023outflows],
where images from protostar _HH211_ have been processed with it.
One of those made the cover of _Nature_ magazine (October 5, 2023, Volume 622 Issue 7981) (Fig. 9).

![Left: image of protostar _HH211_, partially processed with `G'MIC` (cover of Nature, courtesy of [Mark McCaughrean/ESA](https://mastodon.social/@markmccaughrean)). Right: an example of the effect of the `G'MIC` **Banding Denoise** filter on an image of the _IC4553_ galaxy (acquired by the JWST, courtesy of [Judy Schmidt](https://astrodon.social/@spacegeck)).](images/nature.png)

# 5. Conclusions and Perspectives

We presented `G'MIC`, an open-source framework for digital image processing, developed since more than 15 years.
`G'MIC` defines its own scripting language to ease the design and application of image processing pipelines and
to allow the definition of new image processing filters and algorithms.
This is a mature software package that can be used by a wide range of users : experts in image processing,
regular users looking for a tool to quickly retouch or generate images,
or researchers working in fields where image processing can play a role to improve acquired data.
The wide range of available user interfaces makes `G'MIC` a versatile tool for image processing.

# References
