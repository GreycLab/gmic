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
`G'MIC` provides several user interfaces allowing for the manipulation of digital images, adapted to different levels of user expertise,
either from the command line, or as a _C/C`++`_ library, or as a user-friendly graphical plug-in that extends the capabilities of popular digital
image retouching applications, such as _GIMP_, _Krita_, _Photoshop_, _Affinity Photo_ and others.

## Keywords

Image Analysis, Processing and Filtering, Computer Graphics, Scripting Language, User Interfaces, Creative Coding.

# 1. Statement of Need

## 1.1. Context

Intrinsic to `G'MIC`'s design are means to map image processing pipelines to commands, advancing the tool as a self-extending language.
Primal command pipelines may be further assembled into those having wider remits, these suitably named to bespeak their extended purposes and available for succeeding command prototyping.

`G'MIC` is distributed under the CeCILL free software licenses (GPL-compatible). The core language projects several user interfaces to convert,
process or visualize generic *image datasets*. Allied with pipeline toolset, `G'MIC` embodies a highly flexible image model,
ranging from 1D scalar signals to 3D+t sequences of multi-spectral volumetric images, hence including 2D color images.
This makes it a versatile tool for image processing, with a wide range of applications in research, industry and graphic design.

## 1.2. History and Motivation

The `G'MIC` project was started in mid-2008 by research scientists working in the IMAGE team of the _GREYC_ (public research laboratory in France),.
whose area of research is the elaboration of image processing algorithms.

To that end, they first began developing [`CImg`](http://cimg.eu) [@cimg], beginning in 1999 and continuing to the present.
`CImg` is an open-source _C`++`_ library for generic image processing. Here, _generic_ implies a library that addresses structurally diverse imagery: photographs, multi-spectral images (e.g. from satellites), medical images (MRI, X-ray, tomography, etc.) and technical animations, among others.

That said, `CImg` exhibits certain limitations for everyday research work:

1. When one simply wants to apply a predefined algorithm from `CImg` to an image, one needs to write a small, _C`++`_ program. It is only a few lines long, but still it must be compiled before it can be executed. In the context of research work, such mechanics are just so many distractions. The idea of being able to run those algorithms directly from the command line is tempting.

2. Over time, a large number of these small, but purpose-specific, programs has accumulated. They are not broadly useful for integration into the `CImg` library and have become an unruly "collection" of specialized algorithms. By design, they cannot be easily distributed and are difficult to maintain (as opposed to a language having package managers, like Python).

These limitations motivated `G'MIC`'s development, beginning in 2008. Two design objectives came to the fore:

1. Enable _pipelines of image processing algorithms_ that may be directly invoked from the command line, without requiring compilation or linking steps.

2. Gather the implementation of specialized algorithms in a single location, facilitating their evolution, maintenance and distribution.

These objectives, in combination with a desire to write new image processing pipelines and algorithms in the most flexible and concise way possible,
gave rise to the idea of _self-extension_. All these objectives led initially to the development of a specialized scripting language:
the `G'MIC` language, and its associated interpreter, distributed as free software.

## 1.4. Related Software

- **Command-line Interfaces:**

The command line interface `gmic` has been originally inspired by _ImageMagick_ [@imagemagick] and _GraphicsMagick_ [@graphicsmagick], particularly the idea of being able to manipulate digital images from a shell.
The main differences between `G'MIC` and _ImageMagick_/_GraphicsMagick_ are that :

1. The type of images processed is more diverse in `G'MIC`.

2. The possibilities offered by the scripting languages associated with each project, for writing image processing pipelines, are more extensive in `G'MIC` (`G'MIC` language makes it possible to write conditions, loops and multi-threaded pipelines, without having to resort to an external scripting tool, such as `sh` or `bash`).

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

- **`gmic`**, a _command-line_ tool to control the `G'MIC` interpreter from a terminal (Fig. 2).

- **_G'MIC-Qt_** is a _Qt_-based [@qt] graphical interface intended to be used as a _plug-in_ for digital image retouching software, such as _GIMP_, _Krita_, _DigiKam_, _Photoshop_, _Affinity Photo_ and others, or as a _stand-alone_ program (Fig. 3).

- **`G'MIC` Online** is a website where a user can upload a color image and apply one of the _G'MIC-Qt_ filters on it.

- **`libgmic`** and **`libgmic`** are respectively _C`++`_ and _C_ libraries which basically provide a simple _C/C`++`_ _API_ to run a `G'MIC` pipeline on a set of input images.

- **_ZArt_** is a _Qt_-based graphical interface used mainly for demonstration purposes, which applies `G'MIC` filters and effects on streamed webcam images in (almost) real-time.

![The command-line interface `gmic` in action.](images/gmic_cli2.png)

![The _G'MIC-Qt_ plug-in in action.](images/gmic_qt_330.jpg)

## 2.3. Visibility and Community

The `G'MIC` framework has been developed since 2008, mainly in the IMAGE team at the
[_GREYC_ laboratory](https://www.greyc.fr/), a French public research laboratory specialized in computer sciences.
The project web page is [https://gmic.eu](https://gmic.eu).
This website brings together a range of resources, from software download links to documentation and tutorial pages.

The core features of the `G'MIC` interpreter are developed by [David Tschumperlé](https://tschumperle.users.greyc.fr/),
the _G'MIC-Qt_ plug-in by [Sébastien Fourey](https://foureys.users.greyc.fr/), both being permanent researchers at _GREYC_.
The other contributors (for documentation, creation of new filters, or implementation of other user interfaces) can be found on
the [software's forum](https://discuss.pixls.us/c/software/gmic/10), hosted by _Pixls.Us_,
an association which promotes the use of open-source software dedicated to photography and image creation.

The `G'MIC` source code is available on these various github repositories:
[`gmic`](https://github.com/GreycLab/gmic/) (interpreter), [`gmic-qt`](https://github.com/c-koi/gmic-qt/) (plug-in) and
[`gmic-community`](https://github.com/GreycLab/gmic-community/) (external contributions, documentation).

Last but not least, the project provides regular updates on new developments on social networks such as
[Mastodon](https://piaille.fr/@gmic) and [Twitter](https://twitter.com/gmic_eu).

# 3. Examples of Research Work Conducted With `G'MIC`

To illustrate the usefulness of `G'MIC` for research, we list here a few image processing research projects
carried out with `G'MIC`, for algorithm development, prototyping, testing and result generation.
For each, we cite its associated research publication.

- **Patch-Based Image Inpainting**:

`G'MIC` has been used to design and implement an original patch-based image _inpainting_ algorithm in [@buyssens2015exemplar] (Fig. 4).

![Patch-based image inpainting with `G'MIC`. Left: input image. Middle: user-defined mask. Right: inpainting result.](images/inpaint.png)

- **Color LUT Compression**:

We also used `G'MIC` to study the problem of 3D _CLUTs_ compression (Color Look Up Tables), for allowing the efficient storage of
color transformations [@tschumperle2020reconstruction].
`G'MIC` provides more than 1100 _CLUTs_ with approximatively only 4MiB of storage (Fig. 5).

![Principle of the `G'MIC` color _LUT_ compression algorithm. An input _CLUT_ (a) is analyzed and relevant color keypoints are deduced (b) and stored as a small image (c). A perceptual metric is used to ensure that the application of the compressed _CLUT_ on an image is visually similar to the application of the original one.](images/clut_compression2.png)

- **Semi-automatic Colorization of Line Arts**:

Colorizing line art drawings is a problem that illustrators are familiar with, as traditional digital tools (e.g. _Bucket Fill_)
are not always working well, e.g. when lines are anti-aliased or contain gaps in the drawing.
In [@fourey2018fast], we describe a _"Smart coloring"_ algorithm, now implemented in `G'MIC` that analyzes the geometry of the contours
and automatically deduces a reasonable flat-colored layer, from a user-defined set of colored strokes (Fig. 6).

![Color spot extrapolation for automatic lineart colorization.](images/lineart_cat.png)

Note that this colorization algorithm has been also implemented natively in GIMP [@jehan18].

- **Automatic Illumination of Flat-colored Drawings**:

In a similar vein, we have designed an original algorithm that illuminates flat-colorized drawings,
by automatically creating a light and shadow layer from a flat-colored layer [@tschumperle2022automatic] (Fig. 7).

![Principle of the Shape Illumination Algorithm. Left: input image, middle-left: estimated 3D normal map. Right: two examples of different illuminations obtained with the Phong lighting model applied with different parameters.](images/illumination.png)

- **Patch-Based Image Style Transfer**:

Image stylization consists in transforming an input image to give it a pictorial style close to that of a second image (style image).
In 2022, were able to develop a patch-based multi-scale algorithm, with a low algorithmic cost [@samuth2022patch],
now implemented in `G'MIC` (Fig. 8).

![Examples of application of the `G'MIC` style transfer method. An input image (top left) is stylized according to a set of different style images (top row).](images/style_transfer.png)

- **Debanding of Astronomical Images**:

`G'MIC` is known to be used in the astronomy research community, in particular for processing images from the James Webb Space Telescope,
which exhibit band frequency noise, that are well removed with `G'MIC` filter **Banding Denoise**.
`G'MIC` has been cited in [@ray2023outflows], where images from protostar _HH211_ have been processed.
One of those made the cover of _Nature_ of October 2023 (Fig. 9).

![Left: image of protostar _HH211_, partially processed with `G'MIC` (cover of Nature, courtesy of [Mark McCaughrean/ESA](https://mastodon.social/@markmccaughrean)). Right: an example of the effect of the `G'MIC` **Banding Denoise** filter on an image of the _IC4553_ galaxy (acquired by the JWST, courtesy of [Judy Schmidt](https://astrodon.social/@spacegeck)).](images/nature.png)

# References
