/*
 #
 #  File        : gmic_cli.cpp
 #                ( C++ source file )
 #
 #  Description : G'MIC CLI interface - A command-line tool to allow the use
 #                of G'MIC commands from the shell.
 #
 #  Copyright   : David Tschumperlé
 #                ( http://tschumperle.users.greyc.fr/ )
 #
 #  Licenses    : This file is 'dual-licensed', you have to choose one
 #                of the two licenses below to apply.
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
*/

#undef gmic_core
#include <signal.h>
#define cimg_appname "gmic"
#define cimg_library gmic_library
#define CImg gmic_image
#define CImgList gmic_list
#include "CImg.h"
#include "gmic.h"
using namespace gmic_library;

// Fallback function for segfault signals.
#if cimg_OS==1
void gmic_segfault_sigaction(int signal, siginfo_t *si, void *arg) {
  cimg::unused(signal,si,arg);
  cimg::mutex(29);
  std::fprintf(cimg::output(),
               "\n\n%s[gmic] G'MIC encountered a %sfatal error%s%s. "
               "Please submit a bug report, at: %shttps://github.com/GreycLab/gmic/issues%s\n\n",
               cimg::t_red,cimg::t_bold,cimg::t_normal,cimg::t_red,
               cimg::t_bold,cimg::t_normal);
  std::fflush(cimg::output());
  cimg::mutex(29,0);
  std::exit(EXIT_FAILURE);
}
#endif

#if cimg_OS==2
int _CRT_glob = 0; // Disable globbing for msys
#endif

// Main entry
//------------
int main(int argc, char **argv) {

  // Set default output messages stream.
  const bool
    is_debug = cimg_option("-debug",false,0) || cimg_option("debug",false,0),
    is_help = (argc==2 || argc==3) && (!std::strcmp(argv[1],"help") || !std::strcmp(argv[1],"-help") ||
                                       !std::strcmp(argv[1],"h") || !std::strcmp(argv[1],"-h"));
  cimg::output(is_debug?stdout:stderr);

  // Set fallback for segfault signals.
#if cimg_OS==1
  struct sigaction sa;
  std::memset(&sa,0,sizeof(sa));
  sigemptyset(&sa.sa_mask);
  sa.sa_sigaction = gmic_segfault_sigaction;
  sa.sa_flags = SA_SIGINFO;
  sigaction(SIGSEGV,&sa,0);
#endif

  // Init resources folder.
  if (!gmic::init_rc()) {
    std::fprintf(cimg::output(),
                 "\n[gmic] Unable to create resources folder.\n");
    std::fflush(cimg::output());
  }

  // Set special path for curl on Windows
  // (in case the use of libcurl is not enabled).
#if cimg_OS==2
  cimg::curl_path("_gmic\\curl",true);
#endif

  // Declare main G'MIC instance.
  static bool is_abort;
  gmic gmic_instance((char*)0,(char*)0,true,(float*)0,&is_abort,(gmic_pixel_type)0);
  gmic_instance.is_debug = is_debug;
  gmic_instance.set_variable("_host",0,"cli");
  gmic_instance.add_commands("cli_start:");

  // Import update file (from resources directory).
  CImg<char> filename_update, command_updates;
  bool is_invalid_updatefile = false;
  char sep = 0;
  filename_update.assign(1024);
  cimg_snprintf(filename_update,filename_update.width(),"%supdate%u.gmic",
                gmic::path_rc(),gmic_version);
  try { command_updates.load_cimg(filename_update); }
  catch (...) {
    try { command_updates.load_raw(filename_update); }
    catch (...) { }
  }
  if (command_updates) try {
      command_updates.unroll('y');
      command_updates.resize(1,command_updates.height() + 1,1,1,0);
      gmic_instance.add_commands(command_updates);
    } catch (...) { is_invalid_updatefile = true; }
  is_invalid_updatefile|=command_updates && (cimg_sscanf(command_updates," #@gmi%c",&sep)!=1 || sep!='c');
  command_updates.assign();

  // Import user file (in parent of resources directory).
  CImg<char> command_user;
  bool is_invalid_userfile = false;
  const char *const filename_user = gmic::path_user();
  try { command_user.load_raw(filename_user); }
  catch (...) {}
  if (command_user) try {
      command_user.resize(1,command_user.height() + 1,1,1,0);
      gmic_instance.add_commands(command_user,filename_user,is_debug);
    } catch (...) { is_invalid_userfile = true; }
  command_user.assign();

  // Convert 'argv' into G'MIC command line.
  CImgList<char> items;
  if (argc==1) // When no args have been specified
    CImg<char>::string("l[] { cli_noarg onfail }").move_to(items);
  else {
    for (int l = 1; l<argc; ++l) { // Split argv as items
      if (std::strchr(argv[l],' ')) {
        CImg<char>::vector('\"').move_to(items);
        CImg<char>(argv[l],(unsigned int)std::strlen(argv[l])).move_to(items);
        CImg<char>::string("\"").move_to(items);
      } else CImg<char>::string(argv[l]).move_to(items);
      items.back().back()=' ';
    }

    // Determine special mode for running .gmic files as scripts : 'gmic commands.gmic [arguments]'.
    if (argc==2 || argc==3) {
      const char *const ext = cimg::split_filename(argv[1]);
      if (!*ext || !std::strcmp(ext,"gmic")) {
        std::FILE *gmic_file = std::fopen(argv[1],"rb");
        if (gmic_file) {
          bool is_command_file = (bool)*ext;
          if (!*ext) { // In case file has no extension, check it starts with a shebang
            char head[2];
            if (std::fread(head,1,2,gmic_file)==2) {
              std::fseek(gmic_file,0,SEEK_SET);
              if (*head=='#' && head[1]=='!') is_command_file = true;
            }
          }
          if (is_command_file) {
            bool allow_main_ = false;
            gmic gi(0,0,false,0,0,(gmic_pixel_type)0);
            gi.add_commands(gmic_file,argv[1],is_debug,0,0,&allow_main_);
            if (allow_main_ && argc==3) { // Check if command '_main_' has arguments
              const unsigned int hash = (int)gmic::hashcode("_main_",false);
              unsigned int ind = 0;
              if (gmic::search_sorted("_main_",gi.command_names[hash],
                                      gi.command_names[hash].size(),ind)) // Command found
                allow_main_ = (bool)gi.command_has_arguments[hash](ind,0);
            }
            gmic_instance.allow_main_ = allow_main_;
          }
          std::fclose(gmic_file);
        }
      }
    }

    // Determine initial verbosity.
    const char *const s_verbosity = std::getenv("GMIC_VERBOSITY");
    if (!s_verbosity || std::sscanf(s_verbosity,"%d%c",&gmic_instance.verbosity,&sep)!=1)
      gmic_instance.verbosity = gmic_instance.allow_main_?0:is_help?0:
        argc==2 && (!std::strcmp(argv[1],"version") || !std::strcmp(argv[1],"-version"))?0:1;
  }

  // Insert startup command.
  const bool is_first_item_verbose = items.width()>1 &&
    (!std::strncmp("-v ",items[0],3) || !std::strncmp("v ",items[0],2) ||
     !std::strncmp("-verbose ",items[0],9) || !std::strncmp("verbose ",items[0],8));
  items.insert(CImg<char>::string("cli_start , ",false),is_first_item_verbose?2:0);

  if (is_invalid_userfile) { // Display warning message in case of invalid user command file
    CImg<char> tmpstr(1024);
    cimg_snprintf(tmpstr,tmpstr.width(),"warn \"File '\"{/\"%s\"}\"' is not a valid G'MIC command file.\" ",
                  filename_user);
    items.insert(CImg<char>::string(tmpstr.data(),false),is_first_item_verbose?2:0);
  }
  if (is_invalid_updatefile) { // Display warning message in case of invalid user command file
    CImg<char> tmpstr(1024);
    cimg_snprintf(tmpstr,tmpstr.width(),"warn \"File '\"{/\"%s\"}\"' is not a valid G'MIC update file.\" ",
                  filename_update.data());
    items.insert(CImg<char>::string(tmpstr.data(),false),is_first_item_verbose?2:0);
  }

  CImg<char> commands_line(items>'x');
  commands_line.back() = 0;
  items.assign();

  // Launch G'MIC interpreter.
  try {
    CImgList<gmic_pixel_type> images;
    CImgList<char> image_names;
    if (is_help && !cimg::is_file(filename_update.data())) {
      images.insert(gmic::stdlib); CImg<char>::string("stdlib").move_to(image_names);
    }
    gmic_instance.run(commands_line.data(),images,image_names);
  } catch (gmic_exception &e) {
    int error_code = 0;
    bool is_error_code = false;

    const char
      *const it1 = gmic_instance.status?std::strstr(gmic_instance.status,"***"):"",
      *const it2 = it1?std::strstr(it1 + 3,"***"):0;
    if (it2 && std::sscanf(it2,"*** %d%c",&error_code,&sep)!=1) error_code = -1;
    else is_error_code = true;

    if (!is_error_code) {

      // Something went wrong during the pipeline execution.
      if (gmic_instance.verbosity<=0) {
        std::fprintf(cimg::output(),"\n[gmic] %s%s%s%s",
                     cimg::t_red,cimg::t_bold,
                     e.what(),cimg::t_normal);
        std::fflush(cimg::output());
      }
      if (*e.command()) {
        std::fprintf(cimg::output(),"\n[gmic] Command '%s%s%s' has the following description: \n",
                     cimg::t_red,e.command(),cimg::t_normal);
        std::fflush(cimg::output());
        CImgList<gmic_pixel_type> images;
        CImgList<char> image_names;
        images.insert(gmic::stdlib);
        CImg<char>::string("stdlib").move_to(image_names);
        CImg<char> tmp_line(1024);
        cimg_snprintf(tmp_line,tmp_line.width(),
                      "l[] { i raw:\"%s\",uint8 m \"%s\" onfail rm } "
                      "l[] { i raw:\"%s\",uint8 m \"%s\" onfail rm } "
                      "cli_start , "
                      "l[] { $_path_commands foreach { i raw:{n},uint8 rm[0] } onfail rm } "
                      "rv help \"%s\"",
                      filename_update.data(),filename_update.data(),
                      filename_user,filename_user,
                      e.command());
        try {
          gmic(tmp_line,images,image_names);
        } catch (...) { // Fallback in case overloaded version of 'help' crashed
          cimg_snprintf(tmp_line,tmp_line.width(),"help \"%s\"",e.command());
          images.assign().insert(gmic::stdlib);
          image_names.assign();
          gmic(tmp_line,images,image_names);
        }
      } else {
        std::fprintf(cimg::output(),"\n\n");
        std::fflush(cimg::output());
      }
    }
    return error_code;
  }
  return 0;
}
