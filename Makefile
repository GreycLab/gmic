#
#  File        : Makefile
#                ( Makefile for GNU 'make' utility )
#
#  Description : This Makefile exposes the following entries:
#
#                 . 'all':
#                   Equivalent to 'lib','cli','gimp','gmic_qt','libc','zart' (see below).
#
#                 . 'lib':
#                   C++ API for the G'MIC library. Generate files 'src/libgmic.*'.
#
#                 . 'libc' and 'libc_static':
#                   C API for the G'MIC library. Generate files 'src/libcgmic.*'.
#                   'libc_static' tries to embed most libraries as static.
#
#                 . 'cli' and 'cli_shared':
#                    G'MIC command line tool. Generate file 'src/gmic[.exe]'
#                    'cli_shared' generates a binary that is dynamically linked with 'libgmic'.
#
#                 . 'gimp' and 'gimp_shared':
#                   G'MIC-Qt plug-in for GIMP. Generate file 'gmic-qt/gmic_gimp_qt'.
#                   'gimp_shared' generates a binary that is dynamically linked with 'libgmic'.
#
#                 . 'gmic_qt' and 'gmic_qt_shared':
#                   G'MIC-Qt stand-alone application. Generate file 'gmic-qt/gmic_qt'.
#                   'gimp_shared' generates a binary that is dynamically linked with 'libgmic'.
#
#                 . 'zart':
#                   ZArt interface for real-time processing of videos coming from webcams or files.
#                   Generate file 'zart/zart'.
#
#                ( https://gmic.eu )
#
#  Copyright   : David Tschumperlé
#                ( https://tschumperle.users.greyc.fr/ )
#
#  Licenses    : This file is 'dual-licensed', you have to choose one
#      	          of the two licenses below to apply.
#
#                CeCILL-C
#                The CeCILL-C license is close to the GNU LGPL.
#                ( http://cecill.info/licences/Licence_CeCILL-C_V1-en.html )
#
#            or  CeCILL v2.1
#                The CeCILL license is compatible with the GNU GPL.
#                ( http://cecill.info/licences/Licence_CeCILL_V2.1-en.html )
#
#  This software is governed either by the CeCILL or the CeCILL-C license
#  under French law and abiding by the rules of distribution of free software.
#  You can  use, modify and or redistribute the software under the terms of
#  the CeCILL or CeCILL-C licenses as circulated by CEA, CNRS and INRIA
#  at the following URL: "http://cecill.info".
#
#  As a counterpart to the access to the source code and  rights to copy,
#  modify and redistribute granted by the license, users are provided only
#  with a limited warranty  and the software's author,  the holder of the
#  economic rights,  and the successive licensors  have only  limited
#  liability.
#
#  In this respect, the user's attention is drawn to the risks associated
#  with loading,  using,  modifying and/or developing or reproducing the
#  software by the user in light of its specific status of free software,
#  that may mean  that it is complicated to manipulate,  and  that  also
#  therefore means  that it is reserved for developers  and  experienced
#  professionals having in-depth computer knowledge. Users are therefore
#  encouraged to load and test the software's suitability as regards their
#  requirements in conditions enabling the security of their systems and/or
#  data to be ensured and,  more generally, to use and operate it in the
#  same conditions as regards security.
#
#  The fact that you are presently reading this means that you have had
#  knowledge of the CeCILL and CeCILL-C licenses and that you accept its terms.
#

all: all

%:
	cd src && $(MAKE) $*
