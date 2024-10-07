# Contributing to G'MIC

First off, thank you for considering contributing to G'MIC! Your help is highly appreciated. Here are some guidelines to help you get started.

## Table of Contents

- [How Can I Contribute?](#how-can-i-contribute)
  - [Reporting Bugs](#reporting-bugs)
  - [Suggesting Enhancements](#suggesting-enhancements)
  - [Contributing Code](#contributing-code)
- [Style Guides](#style-guides)
  - [Coding Standards](#coding-standards)
  - [Git Commit Messages](#git-commit-messages)
- [Getting General Help](#getting-general-help)

## How Can I Contribute?

### Reporting Bugs

If you find a bug in the G'MIC project, please file a bug report using Github Issues.
Before submitting a report, please check if the issue has already been reported.
Here are the different links you can use, depending on the type of bug to be reported.
Using the right link will optimize the chance of resolution, but if you don't know which one to use, just use the first one:

- [Issues about the G'MIC core functionnalities](https://github.com/dtschump/gmic/issues), i.e. concerning all about the G'MIC language interpreter and the image processing algorithms.
- [Issues about the G'MIC-Qt interface](https://github.com/GreycLab/gmic-qt/issues/), i.e. bugs related to the Qt-based interface of G'MIC.
- [Issues specifically about the 8bf plug-in](https://github.com/0xC0000054/gmic-8bf/) (e.g. installation issues).

If you don't have a github account, you may also consider reporting a bug on the [G'MIC official discussion forum](https://discuss.pixls.us/c/software/gmic/).

#### Bug Report Guidelines
- **Title**: Please provide a clear and descriptive title.
- **Description**: Explain the problem and preferably include steps to reproduce the bug. If we can reproduce the bug on our side, we'll probably already have done most of the work to resolve it.
- **Environment**: Include details about your environment (e.g., OS, G'MIC version).
- **Screenshots**: If applicable, add screenshots or videos to help explain your problem.
- **Logs/Output**: Include any relevant logs or output (e.g. by running G'MIC in the _'Debug'_ mode).

### Suggesting Enhancements

Enhancement suggestions are managed preferably through the [G'MIC official discussion forum](https://discuss.pixls.us/c/software/gmic/).

Please provide a detailed explanation of the enhancement and its potential benefits.

### Contributing Code

Contributing new code to G'MIC is preferably done by submitting new scripts, written in the G'MIC language.
For this, please fill a [Pull Request in the gmic-community repository](https://github.com/GreycLab/gmic-community/pulls).

In case you want to contribute to other aspects of the G'MIC framework (interpreter, native image processing functions, ...), please
submit a PR for the corresponding repository:

- [Improve native (C++) image processing functionalities](https://github.com/GreycLab/CImg/pulls).

Please note that we will not consider including native (C++) contributions in the CImg library if it can be written as a G'MIC script instead.

- [Improve or enhance the G'MIC language interpreter](https://github.com/GreycLab/gmic/pulls).

We're more interested in optimizing the interpreter than in new language features, unless they're essential for executing new G'MIC scripts of major importance (which wouldn't be able to do the same thing without these new features).

- [Improve the G'MIC-Qt interface](https://github.com/GreycLab/gmic-qt/pulls).

Please always open Pull Requests for the `develop` branch of these repositories, when it exists.

## Getting General Help

If you need general help, feel free to ask questions on the [G'MIC discussion forums](https://discuss.pixls.us/c/software/gmic/).
This is where you're most likely to get a precise, detailed answer.

Thank you for your interest in contributing to G'MIC! Your efforts help make the project better for everyone.
