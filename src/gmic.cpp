/*
 #
 #  File        : gmic.cpp
 #                ( C++ source file )
 #
 #  Description : GREYC's Magic for Image Computing - G'MIC core
 #                ( https://gmic.eu )
 #
 #  Copyright   : David Tschumperlé
 #                ( https://tschumperle.users.greyc.fr/ )
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

// Add G'MIC-specific methods to the CImg<T> class of the CImg library.
//----------------------------------------------------------------------
#if defined(cimg_plugin)

template<typename t>
static CImg<T> copy_rounded(const CImg<t>& img) {
  if (!cimg::type<t>::is_float() || cimg::type<T>::is_float()) return img;
  CImg<T> res(img._width,img._height,img._depth,img._spectrum);
  const t *ptrs = img._data;
  cimg_for(res,ptrd,T) *ptrd = (T)cimg::round(*(ptrs++));
  return res;
}

static CImg<T> copy_rounded(const CImg<T>& img) {
  return CImg<T>(img,true);
}

static const char *storage_type(const CImgList<T>& images, const bool allow_bool) {
  T im = cimg::type<T>::max(), iM = cimg::type<T>::min();
  bool is_int = true;
  for (unsigned int l = 0; l<images.size() && is_int; ++l) {
    cimg_for(images[l],p,T) {
      const T val = *p;
      if (!(val==(T)(int)val)) { is_int = false; break; }
      if (val<im) im = val;
      if (val>iM) iM = val;
    }
  }
  if (is_int) {
    if (allow_bool && im==0 && iM==1) return "bool";
    else if (im>=0) {
      if (iM<(1U<<8)) return "uint8";
      else if (iM<(1U<<16)) return "uint16";
      else if (iM<((cimg_uint64)1<<32)) return "uint32";
    } else {
      if (im>=-(1<<7) && iM<(1<<7) && cimg::type<char>::min()<0) return "int8";
      else if (im>=-(1<<15) && iM<(1<<15)) return "int16";
      else if (im>=-((cimg_int64)1<<31) && iM<((cimg_int64)1<<31)) return "int32";
    }
  }
  return pixel_type();
}

static CImg<T> append_CImg3d(const CImgList<T>& images) {
  if (!images) return CImg<T>();
  if (images.size()==1) return +images[0];
  CImg<charT> error_message(1024);
  unsigned int g_nbv = 0, g_nbp = 0;
  ulongT siz = 0;
  cimglist_for(images,l) {
    const CImg<T>& img = images[l];
    if (!img.is_CImg3d(false,error_message))
      throw CImgArgumentException("append_CImg3d(): image [%d] (%u,%u,%u,%u,%p) "
                                  "is not a CImg3d (%s).",
                                  l,img._width,img._height,img._depth,img._spectrum,img._data,
                                  error_message.data());
    siz+=img.size() - 8;
    g_nbv+=cimg::float2uint((float)img[6]);
    g_nbp+=cimg::float2uint((float)img[7]);
  }

  CImg<T> res(1,siz + 8);
  const T **const ptrs = new const T*[images.size()];
  T *ptrd = res._data;
  *(ptrd++) = (T)('C' + 0.5f); *(ptrd++) = (T)('I' + 0.5f); // Create object header
  *(ptrd++) = (T)('m' + 0.5f); *(ptrd++) = (T)('g' + 0.5f);
  *(ptrd++) = (T)('3' + 0.5f); *(ptrd++) = (T)('d' + 0.5f);
  *(ptrd++) = (T)cimg::uint2float(g_nbv);
  *(ptrd++) = (T)cimg::uint2float(g_nbp);
  cimglist_for(images,l) { // Merge object points
    const CImg<T>& img = images[l];
    const unsigned int nbv = cimg::float2uint((float)img[6]);
    std::memcpy(ptrd,img._data + 8,3*nbv*sizeof(T));
    ptrd+=3*nbv;
    ptrs[l] = img._data + 8 + 3*nbv;
  }
  ulongT poff = 0;
  cimglist_for(images,l) { // Merge object primitives
    const unsigned int
      nbv = cimg::float2uint((float)images[l][6]),
      nbp = cimg::float2uint((float)images[l][7]);
    for (unsigned int p = 0; p<nbp; ++p) {
      const unsigned int
        nbi = cimg::float2uint((float)*(ptrs[l]++)),
        _nbi = nbi<5?nbi:nbi==5?2:nbi/3;
      *(ptrd++) = (T)cimg::uint2float(nbi);
      for (unsigned int i = 0; i<_nbi; ++i)
        *(ptrd++) = (T)cimg::uint2float(cimg::float2uint((float)*(ptrs[l]++)) + poff);
      for (unsigned int i = nbi - _nbi; i; --i)
        *(ptrd++) = *(ptrs[l]++);
    }
    poff+=nbv;
  }
  ulongT voff = 0;
  cimglist_for(images,l) { // Merge object colors
    const unsigned int nbc = cimg::float2uint((float)images[l][7]);
    for (unsigned int c = 0; c<nbc; ++c)
      if (*(ptrs[l])==(T)-128) {
        *(ptrd++) = *(ptrs[l]++);
        const unsigned int
          w = (unsigned int)cimg::float2uint((float)*(ptrs[l]++)),
          h = (unsigned int)*(ptrs[l]++),
          s = (unsigned int)*(ptrs[l]++);
        if (!h && !s) { *(ptrd++) = (T)cimg::uint2float((unsigned int)(w + voff)); *(ptrd++) = 0; *(ptrd++) = 0; }
        else {
          *(ptrd++) = (T)w; *(ptrd++) = (T)h; *(ptrd++) = (T)s;
          const ulongT whs = (ulongT)w*h*s;
          std::memcpy(ptrd,ptrs[l],whs*sizeof(T));
          ptrs[l]+=whs; ptrd+=whs;
        }
      } else { *(ptrd++) = *(ptrs[l]++); *(ptrd++) = *(ptrs[l]++); *(ptrd++) = *(ptrs[l]++); }
    voff+=nbc;
  }
  voff = 0;
  cimglist_for(images,l) { // Merge object opacities
    const unsigned int nbo = cimg::float2uint((float)images[l][7]);
    for (unsigned int o = 0; o<nbo; ++o)
      if (*(ptrs[l])==(T)-128) {
        *(ptrd++) = *(ptrs[l]++);
        const unsigned int
          w = (unsigned int)cimg::float2uint((float)*(ptrs[l]++)),
          h = (unsigned int)*(ptrs[l]++),
          s = (unsigned int)*(ptrs[l]++);
        if (!h && !s) { *(ptrd++) = (T)cimg::uint2float((unsigned int)(w + voff)); *(ptrd++) = 0; *(ptrd++) = 0; }
        else {
          *(ptrd++) = (T)w; *(ptrd++) = (T)h; *(ptrd++) = (T)s;
          const ulongT whs = (ulongT)w*h*s;
          std::memcpy(ptrd,ptrs[l],whs*sizeof(T));
          ptrs[l]+=whs; ptrd+=whs;
        }
      } else *(ptrd++) = *(ptrs[l]++);
    voff+=nbo;
  }
  delete[] ptrs;
  return res;
}

CImg<T>& append_string_to(CImg<T>& img, T* &ptrd) const {
  if (!_width) return img;
  if (ptrd + _width>=img.end()) {
    CImg<T> tmp(std::max(8U,2*img._width + _width + 1));
    std::memcpy(tmp,img,img._width*sizeof(T));
    ptrd = tmp._data + (ptrd - img._data);
    tmp.move_to(img);
  }
  std::memcpy(ptrd,_data,_width*sizeof(T));
  ptrd+=_width;
  return img;
}

static CImg<T>& append_string_to(const char c, CImg<T>& img, T* &ptrd) {
  if (ptrd + 1>=img.end()) {
    CImg<T> tmp(std::max(8U,2*img._width + 1));
    std::memcpy(tmp,img,img._width*sizeof(T));
    ptrd = tmp._data + (ptrd - img._data);
    tmp.move_to(img);
  }
  *(ptrd++) = c;
  return img;
}

// Return a copymarked version of an image name.
// This method has no 'in-place' version, at it is always better to call the new instance version.
CImg<T> get_copymark() const {
  if (is_empty() || !*_data) return CImg<T>::string("_c1");

  const char *pe = _data + _width - 1, *ext = cimg::split_filename(_data);
  if (*ext) pe = --ext;
  unsigned int num = 0, fact = 1;
  if (pe>_data + 2) { // Try to find ending number if any
    const char *ppe = pe - 1;
    while (ppe>_data && *ppe>='0' && *ppe<='9') { num+=fact*(*(ppe--) - '0'); fact*=10; }
    if (ppe>_data && ppe!=pe - 1 && *(ppe - 1)=='_' && *ppe=='c' && ppe[1]!='0') {
      pe = ppe - 1;
    }
    else num = 0;
  }
  ++num;
  const unsigned int
    ndigits = (unsigned int)std::max(1.,std::ceil(std::log10(num + 1.))),
    lbase = (unsigned int)(pe - _data),
    lext = _data + _width - ext - 1;
  CImg<T> res(lbase + 2 + ndigits + lext + 1);
  std::memcpy(res,_data,lbase);
  cimg_snprintf(res._data + lbase,res._width - lbase,"_c%u%s",num,ext);
  return res;
}

CImg<T> get_draw_ellipse(const int x, const int y, const float r0, const float r1,
                         const float angle, const T *const col, const float opacity) const {
  return (+*this).draw_ellipse(x,y,r0,r1,angle,col,opacity);
}

CImg<T> get_draw_ellipse(const int x, const int y, const float r0, const float r1,
                         const float angle, const T *const col, const float opacity,
                         const unsigned int pattern) const {
  return (+*this).draw_ellipse(x,y,r0,r1,angle,col,opacity,pattern);
}

CImg<T> get_draw_fill(const int x, const int y, const int z,
                      const T *const col, const float opacity,
                      const float tolerance, const bool is_high_connectivity) const {
  return (+*this).draw_fill(x,y,z,col,opacity,tolerance,is_high_connectivity);
}

template<typename t, typename tc>
CImg<T> get_draw_graph(const CImg<t>& data,
                       const tc *const color, const float opacity,
                       const unsigned int plot_type, const int vertex_type,
                       const double ymin, const double ymax,
                       const unsigned int pattern) const {
  return (+*this).draw_graph(data,color,opacity,plot_type,vertex_type,ymin,ymax,pattern);
}

template<typename t, typename tc>
CImg<T> gmic_draw_graph(const CImg<t>& data,
                        const tc *const color, const float opacity,
                        const unsigned int plot_type, const int vertex_type,
                        const double ymin, const double ymax,
                        const unsigned int pattern) {
  double m = ymin, M = ymax;
  if (ymin==ymax) m = (double)data.max_min(M);
  if (m==M) { --m; ++M; }
  cimg_forC(data,c)
    draw_graph(data.get_shared_channel(c),
               color,opacity,plot_type,vertex_type,m,M,pattern);
  return *this;
}

template<typename t, typename tc>
CImg<T> get_gmic_draw_graph(const CImg<t>& data,
                            const tc *const color, const float opacity,
                            const unsigned int plot_type, const int vertex_type,
                            const double ymin, const double ymax,
                            const unsigned int pattern) const {
  return (+*this).gmic_draw_graph(data,color,opacity,plot_type,vertex_type,ymin,ymax,pattern);
}

CImg<T>& gmic_draw_image(const float x, const float y, const float z, const float c,
                         const char sepx, const char sepy, const char sepz, const char sepc,
                         const CImg<T>& sprite, const CImg<T>& mask, const float opacity,
                         const float max_opacity_mask) {
  const float
    fx = sepx=='~'?x*(width() - sprite.width()):sepx=='%'?x*(width() - 1)/100:x,
    fy = sepy=='~'?y*(height() - sprite.height()):sepy=='%'?y*(height() - 1)/100:y,
    fz = sepz=='~'?y*(depth() - sprite.depth()):sepz=='%'?z*(depth() - 1)/100:z,
    fc = sepc=='~'?c*(spectrum() - sprite.spectrum()):sepc=='%'?c*(spectrum() - 1)/100:c;
  return draw_image((int)cimg::round(fx),(int)cimg::round(fy),
                    (int)cimg::round(fz),(int)cimg::round(fc),
                    sprite,mask,opacity,max_opacity_mask);
}

CImg<T> get_gmic_draw_image(const float x, const float y, const float z, const float c,
                            const char sepx, const char sepy, const char sepz, const char sepc,
                            const CImg<T>& sprite, const CImg<T>& mask, const float opacity,
                            const float max_opacity_mask) const {
  return (+*this).gmic_draw_image(x,y,z,c,sepx,sepy,sepz,sepc,sprite,mask,opacity,max_opacity_mask);
}

CImg<T>& gmic_draw_image(const float x, const float y, const float z, const float c,
                         const char sepx, const char sepy, const char sepz, const char sepc,
                         const CImg<T>& sprite, const float opacity) {
  const float
    fx = sepx=='~'?x*(width() - sprite.width()):sepx=='%'?x*(width() - 1)/100:x,
    fy = sepy=='~'?y*(height() - sprite.height()):sepy=='%'?y*(height() - 1)/100:y,
    fz = sepz=='~'?y*(depth() - sprite.depth()):sepz=='%'?z*(depth() - 1)/100:z,
    fc = sepc=='~'?c*(spectrum() - sprite.spectrum()):sepc=='%'?c*(spectrum() - 1)/100:c;
  return draw_image((int)cimg::round(fx),(int)cimg::round(fy),
                    (int)cimg::round(fz),(int)cimg::round(fc),
                    sprite,opacity);
}

CImg<T> get_gmic_draw_image(const float x, const float y, const float z, const float c,
                            const char sepx, const char sepy, const char sepz, const char sepc,
                            const CImg<T>& sprite, const float opacity) const {
  return (+*this).gmic_draw_image(x,y,z,c,sepx,sepy,sepz,sepc,sprite,opacity);
}

CImg<T> get_draw_line(const int x0, const int y0, const int x1, const int y1, const T *const col,
                      const float opacity, const unsigned int pattern) const {
  return (+*this).draw_line(x0,y0,x1,y1,col,opacity,pattern);
}

CImg<T> get_draw_mandelbrot(const CImg<T>& color_palette, const float opacity,
                            const double z0r, const double z0i, const double z1r, const double z1i,
                            const unsigned int itermax, const bool normalized_iteration,
                            const bool julia_set, const double paramr, const double parami) const {
  return (+*this).draw_mandelbrot(color_palette,opacity,z0r,z0i,z1r,z1i,itermax,
                                  normalized_iteration,julia_set,paramr,parami);
}

template<typename tp, typename tf, typename tc, typename to>
CImg<T> get_draw_object3d(const float x0, const float y0, const float z0,
                          const CImg<tp>& vertices, const CImgList<tf>& primitives,
                          const CImgList<tc>& colors, const CImgList<to>& opacities,
                          const unsigned int render_mode, const bool double_sided,
                          const float focale,
                          const float light_x, const float light_y,const float light_z,
                          const float specular_lightness, const float specular_shininess,
                          const float g_opacity, CImg<floatT>& zbuffer) const {
  return (+*this).draw_object3d(x0,y0,z0,vertices,primitives,colors,opacities,render_mode,
                                double_sided,focale,light_x,light_y,light_z,specular_lightness,
                                specular_shininess,g_opacity,zbuffer);
}

CImg<T> get_draw_plasma(const float alpha, const float beta, const unsigned int scale) const {
  return (+*this).draw_plasma(alpha,beta,scale);
}

CImg<T> get_draw_point(const int x, const int y, const int z, const T *const col,
                       const float opacity) const {
  return (+*this).draw_point(x,y,z,col,opacity);
}

template<typename t>
CImg<T> get_draw_polygon(const CImg<t>& pts, const T *const col, const float opacity) const {
  return (+*this).draw_polygon(pts,col,opacity);
}

template<typename t>
CImg<T> get_draw_polygon(const CImg<t>& pts, const T *const col, const float opacity,
                         const unsigned int pattern) const {
  return (+*this).draw_polygon(pts,col,opacity,pattern);
}

CImg<T>& gmic_autocrop(const CImg<T>& color=CImg<T>::empty()) {
  if (color.width()==1) autocrop(*color);
  else autocrop(color);
  return *this;
}

CImg<T> get_gmic_autocrop(const CImg<T>& color=CImg<T>::empty()) {
  return (+*this).gmic_autocrop(color);
}

CImg<T>& gmic_blur(const float sigma_x, const float sigma_y, const float sigma_z, const float sigma_c,
                   const unsigned int boundary_conditions, const bool is_gaussian) {
  if (is_empty()) return *this;
  if (is_gaussian) {
    if (_width>1) vanvliet(sigma_x,0,'x',boundary_conditions);
    if (_height>1) vanvliet(sigma_y,0,'y',boundary_conditions);
    if (_depth>1) vanvliet(sigma_z,0,'z',boundary_conditions);
    if (_spectrum>1) vanvliet(sigma_c,0,'c',boundary_conditions);
  } else {
    if (_width>1) deriche(sigma_x,0,'x',boundary_conditions);
    if (_height>1) deriche(sigma_y,0,'y',boundary_conditions);
    if (_depth>1) deriche(sigma_z,0,'z',boundary_conditions);
    if (_spectrum>1) deriche(sigma_c,0,'c',boundary_conditions);
  }
  return *this;
}

CImg<Tfloat> get_gmic_blur(const float sigma_x, const float sigma_y, const float sigma_z, const float sigma_c,
                           const unsigned int boundary_conditions, const bool is_gaussian) const {
  return CImg<Tfloat>(*this,false).gmic_blur(sigma_x,sigma_y,sigma_z,sigma_c,boundary_conditions,is_gaussian);
}

CImg<T>& gmic_blur_box(const float sigma_x, const float sigma_y, const float sigma_z, const float sigma_c,
                       const unsigned int order, const unsigned int boundary_conditions,
                       const unsigned int nb_iter) {
  if (is_empty()) return *this;
  if (_width>1) boxfilter(sigma_x,order,'x',boundary_conditions,nb_iter);
  if (_height>1) boxfilter(sigma_y,order,'y',boundary_conditions,nb_iter);
  if (_depth>1) boxfilter(sigma_z,order,'z',boundary_conditions,nb_iter);
  if (_spectrum>1) boxfilter(sigma_c,order,'c',boundary_conditions,nb_iter);
  return *this;
}

CImg<Tfloat> get_gmic_blur_box(const float sigma_x, const float sigma_y, const float sigma_z, const float sigma_c,
                               const unsigned int order, const unsigned int boundary_conditions,
                               const unsigned int nb_iter) const {
  return CImg<Tfloat>(*this,false).gmic_blur_box(sigma_x,sigma_y,sigma_z,sigma_c,order,boundary_conditions,nb_iter);
}

CImg<T>& gmic_blur_box(const float sigma, const unsigned int order, const unsigned int boundary_conditions,
                       const unsigned int nb_iter) {
  const float nsigma = sigma>=0?sigma:-sigma*cimg::max(_width,_height,_depth)/100;
  return gmic_blur_box(nsigma,nsigma,nsigma,0,order,boundary_conditions,nb_iter);
}

CImg<Tfloat> get_gmic_blur_box(const float sigma, const unsigned int order, const unsigned int boundary_conditions,
                               const unsigned int nb_iter) const {
  return CImg<Tfloat>(*this,false).gmic_blur_box(sigma,order,boundary_conditions,nb_iter);
}

CImg<T>& gmic_discard(const char *const axes) {
  for (const char *s = axes; *s; ++s) discard(*s);
  return *this;
}

CImg<T> get_gmic_discard(const char *const axes) const {
  return (+*this).gmic_discard(axes);
}

template<typename t>
CImg<T>& gmic_discard(const CImg<t>& values, const char *const axes) {
  if (is_empty() || !values || !axes || !*axes) return *this;
  for (const char *s = axes; *s; ++s) discard(values,*s);
  return *this;
}

template<typename t>
CImg<T> get_gmic_discard(const CImg<t>& values, const char *const axes) const {
  return (+*this).gmic_discard(values,axes);
}

CImg<T>& gmic_draw_text(const float x, const float y,
                        const char sepx, const char sepy,
                        const char *const text, const T *const col,
                        const int bg, const float opacity, const unsigned int siz,
                        const unsigned int nb_cols) {
  float fx = 0, fy = 0;
  if (is_empty()) {
    const T one[] = { (T)1 };
    fx = sepx=='%' || sepx=='~'?0:x;
    fy = sepy=='%' || sepy=='~'?0:y;
    draw_text((int)cimg::round(fx),(int)cimg::round(fy),"%s",one,0,opacity,siz,text).resize(-100,-100,1,nb_cols);
    cimg_forC(*this,c) if (col[c]!=1) get_shared_channel(c)*=col[c];
    return *this;
  }
  if (sepx=='~' || sepy=='~') {
    const char one[] = { 1 };
    CImg<ucharT> foo;
    foo.draw_text(0,0,"%s",one,0,1,siz,text);
    fx = sepx=='~'?x*(width() - foo.width()):sepx=='%'?x*(width() - 1)/100:x;
    fy = sepy=='~'?y*(height() - foo.height()):sepy=='%'?y*(height() - 1)/100:y;
  } else {
    fx = sepx=='%'?x*(width() - 1)/100:x;
    fy = sepy=='%'?y*(height() - 1)/100:y;
  }
  return draw_text((int)cimg::round(fx),(int)cimg::round(fy),"%s",col,bg,opacity,siz,text);
}

CImg<T> get_gmic_draw_text(const float x, const float y,
                           const char sepx, const char sepy,
                           const char *const text, const T *const col,
                           const int bg, const float opacity, const unsigned int siz,
                           const unsigned int nb_cols) const {
  return (+*this).gmic_draw_text(x,y,sepx,sepy,text,col,bg,opacity,siz,nb_cols);
}

CImg<T>& gmic_draw_text(const float x, const float y,
                        const char sepx, const char sepy,
                        const char *const text, const T *const col,
                        const int bg, const float opacity, const CImgList<T>& font,
                        const unsigned int nb_cols) {
  float fx = 0, fy = 0;
  if (is_empty()) {
    const T one[] = { (T)1 };
    fx = sepx=='%' || sepx=='~'?0:x;
    fy = sepy=='%' || sepy=='~'?0:y;
    draw_text((int)cimg::round(fx),(int)cimg::round(fy),"%s",one,0,opacity,&font,text).resize(-100,-100,1,nb_cols);
    cimg_forC(*this,c) get_shared_channel(c)*=col[c];
    return *this;
  }
  if (sepx=='~' || sepy=='~') {
    const char one[] = { 1 };
    CImg<ucharT> foo;
    foo.draw_text(0,0,"%s",one,0,1,&font,text);
    fx = sepx=='~'?x*(width() - foo.width()):sepx=='%'?x*(width() - 1)/100:x;
    fy = sepy=='~'?y*(height() - foo.height()):sepy=='%'?y*(height() - 1)/100:y;
  } else {
    fx = sepx=='%'?x*(width() - 1)/100:x;
    fy = sepy=='%'?y*(height() - 1)/100:y;
  }
  return draw_text((int)cimg::round(fx),(int)cimg::round(fy),"%s",col,bg,opacity,&font,text);
}

CImg<T> get_gmic_draw_text(const float x, const float y,
                           const char sepx, const char sepy,
                           const char *const text, const T *const col,
                           const int bg, const float opacity, const CImgList<T>& font,
                           const unsigned int nb_cols) const {
  return (+*this).gmic_draw_text(x,y,sepx,sepy,text,col,bg,opacity,font,nb_cols);
}

CImg<T>& gmic_invert_endianness(const char *const stype) {

#define _gmic_invert_endianness(svalue_type,value_type) \
  if (!std::strcmp(stype,svalue_type)) \
    if (pixel_type()==cimg::type<value_type>::string()) invert_endianness(); \
    else CImg<value_type>(*this).invert_endianness().move_to(*this);
  if (!std::strcmp(stype,"bool") ||
      !std::strcmp(stype,"uint8") ||
      !std::strcmp(stype,"int8")) return *this;
  _gmic_invert_endianness("uint16",cimg_uint16)
  else _gmic_invert_endianness("int16",cimg_int16)
    else _gmic_invert_endianness("uint32",cimg_uint32)
      else _gmic_invert_endianness("int32",cimg_int32)
        else _gmic_invert_endianness("uint64",cimg_uint64)
          else _gmic_invert_endianness("int64",cimg_int64)
            else _gmic_invert_endianness("float32",cimg_float32)
              else _gmic_invert_endianness("float64",cimg_float64)
                else invert_endianness();
  return *this;
}

CImg<T> get_gmic_invert_endianness(const char *const stype) const {
  return (+*this).gmic_invert_endianness(stype);
}

CImg<T>& gmic_matchpatch(const CImg<T>& patch_image,
                         const unsigned int patch_width,
                         const unsigned int patch_height,
                         const unsigned int patch_depth,
                         const unsigned int nb_iterations,
                         const unsigned int nb_randoms,
                         const float patch_penalization,
                         const bool is_score,
                         const CImg<T> *const initialization) {
  return get_gmic_matchpatch(patch_image,patch_width,patch_height,patch_depth,
                             nb_iterations,nb_randoms,patch_penalization,is_score,initialization).move_to(*this);
}

CImg<T> get_gmic_matchpatch(const CImg<T>& patch_image,
                            const unsigned int patch_width,
                            const unsigned int patch_height,
                            const unsigned int patch_depth,
                            const unsigned int nb_iterations,
                            const unsigned int nb_randoms,
                            const float patch_penalization,
                            const bool is_score,
                            const CImg<T> *const initialization) const {
  CImg<floatT> score, res;
  res = _matchpatch(patch_image,patch_width,patch_height,patch_depth,
                    nb_iterations,nb_randoms,patch_penalization,
                    initialization?*initialization:CImg<T>::const_empty(),
                    is_score,is_score?score:CImg<floatT>::empty());
  const unsigned int s = res._spectrum;
  if (score) res.resize(-100,-100,-100,s + 1,0).draw_image(0,0,0,s,score);
  return res;
}

const CImg<T>& gmic_print(const char *const title, const bool is_debug,
                          const bool is_valid) const {
  cimg::mutex(29);
  CImg<doubleT> st;
  if (is_valid && !is_empty()) get_stats().move_to(st);
  const ulongT siz = size(), msiz = siz*sizeof(T), siz1 = siz - 1,
    mdisp = msiz<8*1024?0U:msiz<8*1024*1024?1U:2U,
    wh = _width*_height, whd = _width*_height*_depth,
    w1 = _width - 1, wh1 = _width*_height - 1, whd1 = _width*_height*_depth - 1;

  std::fprintf(cimg::output(),"%s%s%s%s:\n  %ssize%s = (%u,%u,%u,%u) [%lu %s of %s%s].\n  %sdata%s = %s",
               cimg::t_magenta,cimg::t_bold,title,cimg::t_normal,
               cimg::t_bold,cimg::t_normal,_width,_height,_depth,_spectrum,
               (unsigned long)(mdisp==0?msiz:(mdisp==1?(msiz>>10):(msiz>>20))),
               mdisp==0?"b":(mdisp==1?"Kio":"Mio"),
               _is_shared?"shared ":"",
               pixel_type(),
               cimg::t_bold,cimg::t_normal,
               is_debug?"":"(");
  if (is_debug) std::fprintf(cimg::output(),"%p = (",(void*)_data);
  if (is_valid) {
    if (is_empty()) std::fprintf(cimg::output(),") [%s].\n",
                                 pixel_type());
    else {
      cimg_foroff(*this,off) {
        std::fprintf(cimg::output(),cimg::type<T>::format_s(),cimg::type<T>::format(_data[off]));
        if (off!=siz1) std::fprintf(cimg::output(),"%s",
                                    off%whd==whd1?" ^ ":
                                    off%wh==wh1?"\\":
                                    off%_width==w1?";":",");
        if (off==11 && siz>24) { off = siz1 - 12; std::fprintf(cimg::output(),"(...),"); }
      }
      std::fprintf(cimg::output(),")%s.\n  %smin%s = %g, %smax%s = %g, %smean%s = %g, "
                   "%sstd%s = %g, %scoords_min%s = (%u,%u,%u,%u), "
                   "%scoords_max%s = (%u,%u,%u,%u).\n",
                   _is_shared?" [shared]":"",
                   cimg::t_bold,cimg::t_normal,st[0],
                   cimg::t_bold,cimg::t_normal,st[1],
                   cimg::t_bold,cimg::t_normal,st[2],
                   cimg::t_bold,cimg::t_normal,std::sqrt(st[3]),
                   cimg::t_bold,cimg::t_normal,(int)st[4],(int)st[5],(int)st[6],(int)st[7],
                   cimg::t_bold,cimg::t_normal,(int)st[8],(int)st[9],(int)st[10],(int)st[11]);
    }
  } else std::fprintf(cimg::output(),"%s%sinvalid pointer%s) [shared %s].\n",
                      cimg::t_red,cimg::t_bold,cimg::t_normal,
                      pixel_type());
  std::fflush(cimg::output());
  cimg::mutex(29,0);
  return *this;
}

CImg<T>& gmic_set(const double value,
                  const int x, const int y, const int z, const int v) {
  (*this).atXYZC(x,y,z,v,(T)0) = (T)value;
  return *this;
}

CImg<T> get_gmic_set(const double value,
                     const int x, const int y, const int z, const int v) const {
  return (+*this).gmic_set(value,x,y,z,v);
}

CImg<T>& gmic_shift(const float delta_x, const float delta_y=0, const float delta_z=0, const float delta_c=0,
                    const unsigned int boundary_conditions=0, const bool interpolation=false) {
  if (is_empty()) return *this;
  const int
    idelta_x = (int)cimg::round(delta_x),
    idelta_y = (int)cimg::round(delta_y),
    idelta_z = (int)cimg::round(delta_z),
    idelta_c = (int)cimg::round(delta_c);
  if (!interpolation ||
      (delta_x==(float)idelta_x && delta_y==(float)idelta_y && delta_z==(float)idelta_z && delta_c==(float)idelta_c))
    return shift(idelta_x,idelta_y,idelta_z,idelta_c,boundary_conditions); // Integer displacement
  return _gmic_shift(delta_x,delta_y,delta_z,delta_c,boundary_conditions).move_to(*this);
}

CImg<T> get_gmic_shift(const float delta_x, const float delta_y=0, const float delta_z=0, const float delta_c=0,
                       const unsigned int boundary_conditions=0, const bool interpolation=false) const {
  if (is_empty()) return CImg<T>::empty();
  const int
    idelta_x = (int)cimg::round(delta_x),
    idelta_y = (int)cimg::round(delta_y),
    idelta_z = (int)cimg::round(delta_z),
    idelta_c = (int)cimg::round(delta_c);
  if (!interpolation ||
      (delta_x==(float)idelta_x && delta_y==(float)idelta_y && delta_z==(float)idelta_z && delta_c==(float)idelta_c))
    return (+*this).shift(idelta_x,idelta_y,idelta_z,idelta_c,boundary_conditions); // Integer displacement
  return _gmic_shift(delta_x,delta_y,delta_z,delta_c,boundary_conditions);
}

CImg<T> _gmic_shift(const float delta_x, const float delta_y=0, const float delta_z=0, const float delta_c=0,
                    const unsigned int boundary_conditions=0) const {
  CImg<T> res(_width,_height,_depth,_spectrum);
  if (delta_c!=0) // 4D shift
    switch (boundary_conditions) {
    case 3 : { // Mirror
      const float w2 = 2.f*width(), h2 = 2.f*height(), d2 = 2.f*depth(), s2 = 2.f*spectrum();
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forXYZC(res,x,y,z,c) {
        const float
          mx = cimg::mod(x - delta_x,w2),
          my = cimg::mod(y - delta_y,h2),
          mz = cimg::mod(z - delta_z,d2),
          mc = cimg::mod(c - delta_c,s2);
        res(x,y,z,c) = _linear_atXYZC(mx<width()?mx:w2 - mx - 1,
                                      my<height()?my:h2 - my - 1,
                                      mz<depth()?mz:d2 - mz - 1,
                                      mc<spectrum()?mc:s2 - mc - 1);
      }
    } break;
    case 2 : // Periodic
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forXYZC(res,x,y,z,c) res(x,y,z,c) = _linear_atXYZC_p(x - delta_x,y - delta_y,z - delta_z,c - delta_c);
      break;
    case 1 : // Neumann
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forXYZC(res,x,y,z,c) res(x,y,z,c) = _linear_atXYZC(x - delta_x,y - delta_y,z - delta_z,c - delta_c);
      break;
    default : // Dirichlet
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forXYZC(res,x,y,z,c) res(x,y,z,c) = linear_atXYZC(x - delta_x,y - delta_y,z - delta_z,c - delta_c,(T)0);
    }
  else if (delta_z!=0) // 3D shift
    switch (boundary_conditions) {
    case 3 : { // Mirror
      const float w2 = 2.f*width(), h2 = 2.f*height(), d2 = 2.f*depth();
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forC(res,c) cimg_forXYZ(res,x,y,z) {
        const float
          mx = cimg::mod(x - delta_x,w2),
          my = cimg::mod(y - delta_y,h2),
          mz = cimg::mod(z - delta_z,d2);
        res(x,y,z,c) = _linear_atXYZ(mx<width()?mx:w2 - mx - 1,
                                     my<height()?my:h2 - my - 1,
                                     mz<depth()?mz:d2 - mz - 1,c);
      }
    } break;
    case 2 : // Periodic
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forC(res,c) cimg_forXYZ(res,x,y,z) res(x,y,z,c) = _linear_atXYZ_p(x - delta_x,y - delta_y,z - delta_z,c);
    break;
    case 1 : // Neumann
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forC(res,c) cimg_forXYZ(res,x,y,z) res(x,y,z,c) = _linear_atXYZ(x - delta_x,y - delta_y,z - delta_z,c);
      break;
    default : // Dirichlet
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forC(res,c) cimg_forXYZ(res,x,y,z) res(x,y,z,c) = linear_atXYZ(x - delta_x,y - delta_y,z - delta_z,c,(T)0);
    }
  else if (delta_y!=0) // 2D shift
    switch (boundary_conditions) {
    case 3 : { // Mirror
      const float w2 = 2.f*width(), h2 = 2.f*height();
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forZC(res,z,c) cimg_forXY(res,x,y) {
        const float
          mx = cimg::mod(x - delta_x,w2),
          my = cimg::mod(y - delta_y,h2);
        res(x,y,z,c) = _linear_atXY(mx<width()?mx:w2 - mx - 1,
                                    my<height()?my:h2 - my - 1,z,c);
      }
    } break;
    case 2 : // Periodic
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forZC(res,z,c) cimg_forXY(res,x,y) res(x,y,z,c) = _linear_atXY_p(x - delta_x,y - delta_y,z,c);
      break;
    case 1 : // Neumann
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forZC(res,z,c) cimg_forXY(res,x,y) res(x,y,z,c) = _linear_atXY(x - delta_x,y - delta_y,z,c);
      break;
    default : // Dirichlet
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forZC(res,z,c) cimg_forXY(res,x,y) res(x,y,z,c) = linear_atXY(x - delta_x,y - delta_y,z,c,(T)0);
    }
  else // 1D shift
    switch (boundary_conditions) {
    case 3 : { // Mirror
      const float w2 = 2.f*width();
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forYZC(res,y,z,c) cimg_forX(res,x) {
        const float mx = cimg::mod(x - delta_x,w2);
        res(x,y,z,c) = _linear_atX(mx<width()?mx:w2 - mx - 1,y,z,c);
      }
    } break;
    case 2 : // Periodic
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forYZC(res,y,z,c) cimg_forX(res,x) res(x,y,z,c) = _linear_atX_p(x - delta_x,y,z,c);
      break;
    case 1 : // Neumann
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forYZC(res,y,z,c) cimg_forX(res,x) res(x,y,z,c) = _linear_atX(x - delta_x,y,z,c);
      break;
    default : // Dirichlet
      cimg_pragma_openmp(parallel for cimg_openmp_collapse(3) cimg_openmp_if_size(res.size(),4096))
      cimg_forYZC(res,y,z,c) cimg_forX(res,x) res(x,y,z,c) = linear_atX(x - delta_x,y,z,c,(T)0);
    }
  return res;
}

template<typename t>
const CImg<T>& gmic_symmetric_eigen(CImg<t>& val, CImg<t>& vec) const {
  if (spectrum()!=3 && spectrum()!=6) return symmetric_eigen(val,vec);
  val.assign(width(),height(),depth(),spectrum()==3?2:3);
  vec.assign(width(),height(),depth(),spectrum()==3?2:6);
  CImg<t> _val, _vec;
  cimg_forXYZ(*this,x,y,z) {
    get_tensor_at(x,y,z).symmetric_eigen(_val,_vec);
    val.set_vector_at(_val,x,y,z);
    if (spectrum()==3) {
      vec(x,y,z,0) = _vec(0,0);
      vec(x,y,z,1) = _vec(0,1);
    } else {
      vec(x,y,z,0) = _vec(0,0);
      vec(x,y,z,1) = _vec(0,1);
      vec(x,y,z,2) = _vec(0,2);

      vec(x,y,z,3) = _vec(1,0);
      vec(x,y,z,4) = _vec(1,1);
      vec(x,y,z,5) = _vec(1,2);
    }
  }
  return *this;
}

template<typename t>
CImg<T>& inpaint(const CImg<t>& mask, const unsigned int method) {
  if (!is_sameXYZ(mask))
    throw CImgArgumentException("CImg<%s>::inpaint(): Invalid mask (%u,%u,%u,%u,%p) for "
                                "instance image (%u,%u,%u,%u,%p).",
                                pixel_type(),mask._width,mask._height,mask._depth,
                                mask._spectrum,mask._data,
                                _width,_height,_depth,_spectrum,_data);
  CImg<t> _mask(mask,false), _nmask(mask,false);
  bool is_pixel = false;

  do {
    is_pixel = false;

    if (depth()==1) { // 2D image
      CImg_3x3(M,t);
      CImg_3x3(I,T);

      switch (method) {
      case 0: // Average 2D (low-connectivity)
        cimg_for3x3(_mask,x,y,0,0,M,t) if (Mcc && (!Mcp || !Mpc || !Mnc || !Mcn)) {
          is_pixel = true;
          const unsigned int wcp = Mcp?0U:1U, wpc = Mpc?0U:1U, wnc = Mnc?0U:1U, wcn = Mcn?0U:1U,
            sumw = wcp + wpc + wnc + wcn;
          cimg_forC(*this,k) {
            cimg_get3x3(*this,x,y,0,k,I,T);
            (*this)(x,y,k) = (T)((wcp*Icp + wpc*Ipc + wnc*Inc + wcn*Icn)/(float)sumw);
          }
          _nmask(x,y) = 0;
        }
        break;

      case 1: // Average 2D (high-connectivity)
        cimg_for3x3(_mask,x,y,0,0,M,t) if (Mcc && (!Mpp || !Mcp || !Mnp || !Mpc || !Mnc || !Mpn || !Mcn || !Mnn)) {
          is_pixel = true;
          const unsigned int
            wpp = Mpp?0U:1U, wcp = Mcp?0U:2U, wnp = Mnp?0U:1U,
            wpc = Mpc?0U:2U, wnc = Mnc?0U:2U,
            wpn = Mpn?0U:1U, wcn = Mcn?0U:2U, wnn = Mnn?0U:1U,
            sumw = wpp + wcp + wnp + wpc + wnc + wpn + wcn + wnn;
          cimg_forC(*this,k) {
            cimg_get3x3(*this,x,y,0,k,I,T);
            (*this)(x,y,k) = (T)((wpp*Ipp + wcp*Icp + wnp*Inp + wpc*Ipc +
                                  wnc*Inc + wpn*Ipn + wcn*Icn + wnn*Inn)/(float)sumw);
          }
          _nmask(x,y) = 0;
        }
        break;

      case 2: { // Median 2D (low-connectivity)
        T J[4];
        cimg_for3x3(_mask,x,y,0,0,M,t)
          if (Mcc && (!Mcp || !Mpc || !Mnc || !Mcn)) {
            is_pixel = true;
            cimg_forC(*this,k) {
              cimg_get3x3(*this,x,y,0,k,I,T);
              unsigned int ind = 0;
              if (!Mcp) J[ind++] = Icp;
              if (!Mpc) J[ind++] = Ipc;
              if (!Mnc) J[ind++] = Inc;
              if (!Mcn) J[ind++] = Icn;
              (*this)(x,y,k) = CImg<T>(J,ind,1,1,1,true).kth_smallest(ind>>1);
            }
            _nmask(x,y) = 0;
          }
      } break;

      default: // Median 2D (high-connectivity)
        T J[8];
        cimg_for3x3(_mask,x,y,0,0,M,t)
          if (Mcc && (!Mpp || !Mcp || !Mnp || !Mpc || !Mnc || !Mpn || !Mcn || !Mnn)) {
            is_pixel = true;
            cimg_forC(*this,k) {
              cimg_get3x3(*this,x,y,0,k,I,T);
              unsigned int ind = 0;
              if (!Mpp) J[ind++] = Ipp;
              if (!Mcp) J[ind++] = Icp;
              if (!Mnp) J[ind++] = Inp;
              if (!Mpc) J[ind++] = Ipc;
              if (!Mnc) J[ind++] = Inc;
              if (!Mpn) J[ind++] = Ipn;
              if (!Mcn) J[ind++] = Icn;
              if (!Mnn) J[ind++] = Inn;
              (*this)(x,y,k) = CImg<T>(J,ind,1,1,1,true).kth_smallest(ind>>1);
            }
            _nmask(x,y) = 0;
          }
      }

    } else { // 3D image
      CImg_3x3x3(M,t);
      CImg_3x3x3(I,T);

      switch (method) {
      case 0: // Average 3D (low-connectivity)
        cimg_for3x3x3(_mask,x,y,z,0,M,t)
          if (Mccc && (!Mccp || !Mcpc || !Mpcc || !Mncc || !Mcnc || !Mccn)) {
            is_pixel = true;
            const unsigned int
              wccp = Mccp?0U:1U, wcpc = Mcpc?0U:1U, wpcc = Mpcc?0U:1U,
              wncc = Mncc?0U:1U, wcnc = Mcnc?0U:1U, wccn = Mccn?0U:1U,
              sumw = wcpc + wpcc + wccp + wncc + wcnc + wccn;
            cimg_forC(*this,k) {
              cimg_get3x3x3(*this,x,y,z,k,I,T);
              (*this)(x,y,z,k) = (T)((wccp*Iccp + wcpc*Icpc + wpcc*Ipcc +
                                      wncc*Incc + wcnc*Icnc + wccn*Iccn)/(float)sumw);
            }
            _nmask(x,y,z) = 0;
          }
        break;

      case 1: // Average 3D (high-connectivity)
        cimg_for3x3x3(_mask,x,y,z,0,M,t)
          if (Mccc && (!Mppp || !Mcpp || !Mnpp || !Mpcp || !Mccp || !Mncp || !Mpnp || !Mcnp ||
                       !Mnnp || !Mppc || !Mcpc || !Mnpc || !Mpcc || !Mncc || !Mpnc || !Mcnc ||
                       !Mnnc || !Mppn || !Mcpn || !Mnpn || !Mpcn || !Mccn || !Mncn || !Mpnn ||
                       !Mcnn || !Mnnn)) {
            is_pixel = true;
            const unsigned int
              wppp = Mppp?0U:1U, wcpp = Mcpp?0U:2U, wnpp = Mnpp?0U:1U,
              wpcp = Mpcp?0U:2U, wccp = Mccp?0U:4U, wncp = Mncp?0U:2U,
              wpnp = Mpnp?0U:1U, wcnp = Mcnp?0U:2U, wnnp = Mnnp?0U:1U,
              wppc = Mppc?0U:2U, wcpc = Mcpc?0U:4U, wnpc = Mnpc?0U:2U,
              wpcc = Mpcc?0U:4U, wncc = Mncc?0U:4U,
              wpnc = Mpnc?0U:2U, wcnc = Mcnc?0U:4U, wnnc = Mnnc?0U:2U,
              wppn = Mppn?0U:1U, wcpn = Mcpn?0U:2U, wnpn = Mnpn?0U:1U,
              wpcn = Mpcn?0U:2U, wccn = Mccn?0U:4U, wncn = Mncn?0U:2U,
              wpnn = Mpnn?0U:1U, wcnn = Mcnn?0U:2U, wnnn = Mnnn?0U:1U,
              sumw = wppp + wcpp + wnpp + wpcp + wccp + wncp + wpnp + wcnp + wnnp +
              wppc + wcpc + wnpc + wpcc + wncc + wpnc + wcnc + wnnc +
              wppn + wcpn + wnpn + wpcn + wccn + wncn + wpnn + wcnn + wnnn;
            cimg_forC(*this,k) {
              cimg_get3x3x3(*this,x,y,z,k,I,T);
              (*this)(x,y,z,k) = (T)((wppp*Ippp + wcpp*Icpp + wnpp*Inpp +
                                      wpcp*Ipcp + wccp*Iccp + wncp*Incp +
                                      wpnp*Ipnp + wcnp*Icnp + wnnp*Innp +
                                      wppc*Ippc + wcpc*Icpc + wnpc*Inpc +
                                      wpcc*Ipcc + wncc*Incc +
                                      wpnc*Ipnc + wcnc*Icnc + wnnc*Innc +
                                      wppn*Ippn + wcpn*Icpn + wnpn*Inpn +
                                      wpcn*Ipcn + wccn*Iccn + wncn*Incn +
                                      wpnn*Ipnn + wcnn*Icnn + wnnn*Innn)/(float)sumw);
            }
            _nmask(x,y,z) = 0;
          }
        break;

      case 2: { // Median 3D (low-connectivity)
        T J[6];
        cimg_for3x3x3(_mask,x,y,z,0,M,t)
          if (Mccc && (!Mccp || !Mcpc || !Mpcc || !Mncc || !Mcnc || !Mccn)) {
            is_pixel = true;
            cimg_forC(*this,k) {
              cimg_get3x3x3(*this,x,y,z,k,I,T);
              unsigned int ind = 0;
              if (!Mccp) J[ind++] = Iccp;
              if (!Mcpc) J[ind++] = Icpc;
              if (!Mpcc) J[ind++] = Ipcc;
              if (!Mncc) J[ind++] = Incc;
              if (!Mcnc) J[ind++] = Icnc;
              if (!Mccn) J[ind++] = Iccn;
              (*this)(x,y,z,k) = CImg<T>(J,ind,1,1,1,true).kth_smallest(ind>>1);
            }
            _nmask(x,y,z) = 0;
          }
      } break;

      default: { // Median 3D (high-connectivity)
        T J[26];
        cimg_for3x3x3(_mask,x,y,z,0,M,t)
          if (Mccc && (!Mppp || !Mcpp || !Mnpp || !Mpcp || !Mccp || !Mncp || !Mpnp || !Mcnp ||
                       !Mnnp || !Mppc || !Mcpc || !Mnpc || !Mpcc || !Mncc || !Mpnc || !Mcnc ||
                       !Mnnc || !Mppn || !Mcpn || !Mnpn || !Mpcn || !Mccn || !Mncn || !Mpnn ||
                       !Mcnn || !Mnnn)) {
            is_pixel = true;
            cimg_forC(*this,k) {
              cimg_get3x3x3(*this,x,y,z,k,I,T);
              unsigned int ind = 0;
              if (!Mppp) J[ind++] = Ippp;
              if (!Mcpp) J[ind++] = Icpp;
              if (!Mnpp) J[ind++] = Inpp;
              if (!Mpcp) J[ind++] = Ipcp;
              if (!Mccp) J[ind++] = Iccp;
              if (!Mncp) J[ind++] = Incp;
              if (!Mpnp) J[ind++] = Ipnp;
              if (!Mcnp) J[ind++] = Icnp;
              if (!Mnnp) J[ind++] = Innp;
              if (!Mppc) J[ind++] = Ippc;
              if (!Mcpc) J[ind++] = Icpc;
              if (!Mnpc) J[ind++] = Inpc;
              if (!Mpcc) J[ind++] = Ipcc;
              if (!Mncc) J[ind++] = Incc;
              if (!Mpnc) J[ind++] = Ipnc;
              if (!Mcnc) J[ind++] = Icnc;
              if (!Mnnc) J[ind++] = Innc;
              if (!Mppn) J[ind++] = Ippn;
              if (!Mcpn) J[ind++] = Icpn;
              if (!Mnpn) J[ind++] = Inpn;
              if (!Mpcn) J[ind++] = Ipcn;
              if (!Mccn) J[ind++] = Iccn;
              if (!Mncn) J[ind++] = Incn;
              if (!Mpnn) J[ind++] = Ipnn;
              if (!Mcnn) J[ind++] = Icnn;
              if (!Mnnn) J[ind++] = Innn;
              (*this)(x,y,z,k) = CImg<T>(J,ind,1,1,1,true).kth_smallest(ind>>1);
            }
            _nmask(x,y,z) = 0;
          }
      } break;
      }
    }

    _mask = _nmask;
  } while (is_pixel);
  return *this;
}

template<typename t>
CImg<T> get_inpaint(const CImg<t>& mask, const unsigned int method) const {
  return (+*this).inpaint(mask,method);
}

template<typename t>
CImg<T>& inpaint_patch(const CImg<t>& mask, const unsigned int patch_size=11,
                       const unsigned int lookup_size=22, const float lookup_factor=1,
                       const int lookup_increment=1,
                       const unsigned int blend_size=0, const float blend_threshold=0.5f,
                       const float blend_decay=0.02f, const unsigned int blend_scales=10,
                       const bool is_blend_outer=false) {
  if (depth()>1)
    throw CImgInstanceException(_cimg_instance
                                "inpaint_patch(): Instance image is volumetric (should be 2D).",
                                cimg_instance);
  if (!is_sameXYZ(mask))
    throw CImgArgumentException(_cimg_instance
                                "inpaint_patch() : Sizes of instance image and specified mask "
                                "(%u,%u,%u,%u) do not match.",
                                cimg_instance,
                                mask._width,mask._height,mask._depth,mask._spectrum);
  if (!patch_size)
    throw CImgArgumentException(_cimg_instance
                                "inpaint_patch() : Specified patch size is 0, must be strictly "
                                "positive.",
                                cimg_instance);
  if (!lookup_size)
    throw CImgArgumentException(_cimg_instance
                                "inpaint_patch() : Specified lookup size is 0, must be strictly "
                                "positive.",
                                cimg_instance);
  if (lookup_factor<0)
    throw CImgArgumentException(_cimg_instance
                                "inpaint_patch() : Specified lookup factor %g is negative, must be "
                                "positive.",
                                cimg_instance,
                                lookup_factor);
  if (!lookup_increment)
    throw CImgArgumentException(_cimg_instance
                                "inpaint_patch() : Specified lookup increment is 0, must be "
                                "strictly positive.",
                                cimg_instance);
  if (blend_decay<0)
    throw CImgArgumentException(_cimg_instance
                                "inpaint_patch() : Specified blend decay %g is negative, must be "
                                "positive.",
                                cimg_instance,
                                blend_decay);

  // Find (dilated by 2) bounding box for the inpainting mask.
  unsigned int xm0 = _width, ym0 = _height, xm1 = 0, ym1 = 0;
  bool is_mask_found = false;
  cimg_forXY(mask,x,y) if (mask(x,y)) {
    is_mask_found = true;
    if (x<(int)xm0) xm0 = (unsigned int)x;
    if (x>(int)xm1) xm1 = (unsigned int)x;
    if (y<(int)ym0) ym0 = (unsigned int)y;
    if (y>(int)ym1) ym1 = (unsigned int)y;
  }
  if (!is_mask_found) return *this;
  xm0 = xm0>2?xm0 - 2:0;
  ym0 = ym0>2?ym0 - 2:0;
  xm1 = xm1<_width - 3?xm1 + 2:_width - 1;
  ym1 = ym1<_height - 3?ym1 + 2:_height - 1;
  int ox = (int)xm0, oy = (int)ym0;
  unsigned int dx = xm1 - xm0 + 1U, dy = ym1 - ym0 + 1U;

  // Construct normalized version of the mask.
  CImg<ucharT> nmask(dx,dy);
  {
    unsigned char *ptrM = nmask.data();
    cimg_for_inXY(mask,xm0,ym0,xm1,ym1,x,y) *(ptrM++) = mask(x,y)?0U:1U;
  }
  xm0 = ym0 = 0; xm1 = dx - 1; ym1 = dy - 1;

  // Start patch filling algorithm.
  const int p2 = (int)patch_size/2, p1 = (int)patch_size - p2 - 1;
  const unsigned int patch_size2 = patch_size*patch_size;
  unsigned int _lookup_size = lookup_size, nb_lookups = 0, nb_fails = 0, nb_saved_patches = 0;
  bool is_strict_search = true;
  const float one = 1;

  CImg<floatT> confidences(nmask), priorities(dx,dy,1,2,-1), pC;
  CImg<unsigned int> saved_patches(4,256), is_visited(width(),height(),1,1,0);
  CImg<ucharT> pM, pN;  // Pre-declare patch variables (avoid iterative memory alloc/dealloc)
  CImg<T> pP, pbest;
  CImg<floatT> weights(patch_size,patch_size,1,1,0);
  weights.draw_gaussian((float)p1,(float)p1,patch_size/15.f,&one)/=patch_size2;
  unsigned int target_index = 0;

  while (true) {

    // Extract mask border points and compute priorities to find target point.
    unsigned int nb_border_points = 0;
    float target_confidence = -1, target_priority = -1;
    int target_x = -1, target_y = -1;
    CImg_5x5(M,unsigned char);

    cimg_for_in5x5(nmask,xm0,ym0,xm1,ym1,x,y,0,0,M,unsigned char)
      if (!Mcc && (Mcp || Mcn || Mpc || Mnc)) { // Found mask border point

        float confidence_term = -1, data_term = -1;
        if (priorities(x,y)>=0) { // If priority has already been computed
          confidence_term = priorities(x,y,0);
          data_term = priorities(x,y,1);
        } else { // If priority must be computed/updated

          // Compute smoothed normal vector.
          const float
            // N = smoothed 3x3 neighborhood of M.
            Npc = (4.f*Mpc + 2.f*Mbc + 2.f*Mcc + 2.f*Mpp + 2.f*Mpn + Mbp + Mbn + Mcp + Mcn)/16,
            Nnc = (4.f*Mnc + 2.f*Mac + 2.f*Mcc + 2.f*Mnp + 2.f*Mnn + Map + Man + Mcp + Mcn)/16,
            Ncp = (4.f*Mcp + 2.f*Mcb + 2.f*Mcc + 2.f*Mpp + 2.f*Mnp + Mpb + Mnb + Mpc + Mnc)/16,
            Ncn = (4.f*Mcn + 2.f*Mca + 2.f*Mcc + 2.f*Mpn + 2.f*Mnn + Mpa + Mna + Mpc + Mnc)/16,
            _nx = 0.5f*(Nnc - Npc),
            _ny = 0.5f*(Ncn - Ncp),
            nn = std::sqrt(1e-8f + _nx*_nx + _ny*_ny),
            nx = _nx/nn,
            ny = _ny/nn;

          // Compute confidence term.
          nmask._inpaint_patch_crop(x - p1,y - p1,x + p2,y + p2,1).move_to(pM);
          confidences._inpaint_patch_crop(x - p1,y - p1,x + p2,y + p2,1).move_to(pC);
          confidence_term = 0;
          const unsigned char *ptrM = pM.data();
          cimg_for(pC,ptrC,float) confidence_term+=*ptrC**(ptrM++);
          confidence_term/=patch_size2;
          priorities(x,y,0) = confidence_term;

          // Compute data term.
          _inpaint_patch_crop(ox + x - p1,oy + y - p1,ox + x + p2,oy + y + p2,2).move_to(pP);
          float mean_ix2 = 0, mean_ixiy = 0, mean_iy2 = 0;

          CImg_3x3(I,T);
          CImg_3x3(_M, unsigned char);
          cimg_forC(pP,c) cimg_for3x3(pP,p,q,0,c,I,T) {
            // Compute weight-mean of structure tensor inside patch.
            cimg_get3x3(pM,p,q,0,0,_M,unsigned char);
            const float
              ixf = (float)(_Mnc*_Mcc*(Inc - Icc)),
              iyf = (float)(_Mcn*_Mcc*(Icn - Icc)),
              ixb = (float)(_Mcc*_Mpc*(Icc - Ipc)),
              iyb = (float)(_Mcc*_Mcp*(Icc - Icp)),
              ix = cimg::abs(ixf)>cimg::abs(ixb)?ixf:ixb,
              iy = cimg::abs(iyf)>cimg::abs(iyb)?iyf:iyb,
              w = weights(p,q);
            mean_ix2 += w*ix*ix;
            mean_ixiy += w*ix*iy;
            mean_iy2 += w*iy*iy;
          }
          const float // Compute tensor-directed data term
            ux = mean_ix2*(-ny) + mean_ixiy*nx,
            uy = mean_ixiy*(-ny) + mean_iy2*nx;
          data_term = std::sqrt(ux*ux + uy*uy);
          priorities(x,y,1) = data_term;
        }
        const float priority = confidence_term*data_term;
        if (priority>target_priority) {
          target_priority = priority; target_confidence = confidence_term;
          target_x = ox + x; target_y = oy + y;
        }
        ++nb_border_points;
      }
    if (!nb_border_points) break; // No more mask border points to inpaint!

    // Locate already reconstructed neighbors (if any), to get good origins for patch lookup.
    CImg<unsigned int> lookup_candidates(2,256);
    unsigned int nb_lookup_candidates = 0, *ptr_lookup_candidates = lookup_candidates.data();
    {
      const unsigned int *ptr_saved_patches = saved_patches.data();
      const int
        x0 = target_x - (int)patch_size, y0 = target_y - (int)patch_size,
        x1 = target_x + (int)patch_size, y1 = target_y + (int)patch_size;
      for (unsigned int k = 0; k<nb_saved_patches; ++k) {
        const unsigned int
          src_x = *(ptr_saved_patches++), src_y = *(ptr_saved_patches++),
          dest_x = *(ptr_saved_patches++), dest_y = *(ptr_saved_patches++);
        if ((int)dest_x>=x0 && (int)dest_y>=y0 && (int)dest_x<=x1 && (int)dest_y<=y1) {
          const int off_x = target_x - (int)dest_x, off_y = target_y - (int)dest_y;
          *(ptr_lookup_candidates++) = src_x + off_x;
          *(ptr_lookup_candidates++) = src_y + off_y;
          if (++nb_lookup_candidates>=lookup_candidates._height) {
            lookup_candidates.resize(2,-200,1,1,0);
            ptr_lookup_candidates = lookup_candidates.data(0,nb_lookup_candidates);
          }
        }
      }
    }
    // Add also target point as a center for the patch lookup.
    if (++nb_lookup_candidates>=lookup_candidates._height) {
      lookup_candidates.resize(2,-200,1,1,0);
      ptr_lookup_candidates = lookup_candidates.data(0,nb_lookup_candidates);
    }
    *(ptr_lookup_candidates++) = (unsigned int)target_x;
    *(ptr_lookup_candidates++) = (unsigned int)target_y;

    // Divide size of lookup regions if several lookup sources have been detected.
    unsigned int final_lookup_size = _lookup_size;
    if (nb_lookup_candidates>1) {
      const unsigned int
        _final_lookup_size = std::max(5U,(unsigned int)cimg::round(_lookup_size*lookup_factor/
                                                                    std::sqrt((float)nb_lookup_candidates),1,1));
      final_lookup_size = _final_lookup_size + 1 - (_final_lookup_size%2);
    }
    const int l2 = (int)final_lookup_size/2, l1 = (int)final_lookup_size - l2 - 1;

    // Find best patch candidate to fill target point.
    _inpaint_patch_crop(target_x - p1,target_y - p1,target_x + p2,target_y + p2,0).move_to(pP);
    nmask._inpaint_patch_crop(target_x - ox - p1,target_y - oy - p1,target_x - ox + p2,target_y - oy + p2,0).
      move_to(pM);
    ++target_index;
    const unsigned int
      _lookup_increment = (unsigned int)(lookup_increment>0?lookup_increment:
                                         nb_lookup_candidates>1?1:-lookup_increment);
    float best_ssd = cimg::type<float>::max();
    int best_x = -1, best_y = -1;
    for (unsigned int C = 0; C<nb_lookup_candidates; ++C) {
      const int
        xl = (int)lookup_candidates(0U,C),
        yl = (int)lookup_candidates(1U,C),
        xl0 = std::max(p1,xl - l1), yl0 = std::max(p1,yl - l1),
        xl1 = std::min(width() - 1 - p2,xl + l2), yl1 = std::min(height() - 1 - p2,yl + l2);
      for (int y = yl0; y<=yl1; y+=_lookup_increment)
        for (int x = xl0; x<=xl1; x+=_lookup_increment) if (is_visited(x,y)!=target_index) {
            if (is_strict_search) mask._inpaint_patch_crop(x - p1,y - p1,x + p2,y + p2,1).move_to(pN);
            else nmask._inpaint_patch_crop(x - ox - p1,y - oy - p1,x - ox + p2,y - oy + p2,0).move_to(pN);
            if ((is_strict_search && pN.sum()==0) || (!is_strict_search && pN.sum()==patch_size2)) {
              _inpaint_patch_crop(x - p1,y - p1,x + p2,y + p2,0).move_to(pC);
              float ssd = 0;
              const T *_pP = pP._data;
              const float *_pC = pC._data;
              cimg_for(pM,_pM,unsigned char) { if (*_pM) {
                  cimg_forC(pC,c) {
                    ssd+=cimg::sqr((Tfloat)*_pC - (Tfloat)*_pP); _pC+=patch_size2; _pP+=patch_size2;
                  }
                  if (ssd>=best_ssd) break;
                  _pC-=pC._spectrum*patch_size2;
                  _pP-=pC._spectrum*patch_size2;
                }
                ++_pC; ++_pP;
              }
              if (ssd<best_ssd) { best_ssd = ssd; best_x = x; best_y = y; }
            }
            is_visited(x,y) = target_index;
          }
    }

    if (best_x<0) { // If no best patch found
      priorities(target_x - ox,target_y - oy,0)/=10; // Reduce its priority (lower data_term)
      if (++nb_fails>=4) { // If too much consecutive fails :
        nb_fails = 0;
        _lookup_size+=_lookup_size/2; // Try to expand the lookup size
        if (++nb_lookups>=3) {
          if (is_strict_search) { // If still fails, switch to non-strict search mode
            is_strict_search = false;
            _lookup_size = lookup_size;
            nb_lookups = 0;
          }
          else return *this; // Pathological case, probably a weird mask
        }
      }
    } else { // Best patch found -> reconstruct missing part on the target patch
      _lookup_size = lookup_size;
      nb_lookups = nb_fails = 0;
      _inpaint_patch_crop(best_x - p1,best_y - p1,best_x + p2,best_y + p2,0).move_to(pbest);
      nmask._inpaint_patch_crop(target_x - ox - p1,target_y - oy - p1,target_x - ox + p2,target_y - oy + p2,1).
        move_to(pM);
      cimg_for(pM,ptr,unsigned char) *ptr = (unsigned char)(1 - *ptr);
      draw_image(target_x - p1,target_y - p1,pbest,pM,1,1);
      confidences.draw_image(target_x - ox - p1,target_y - oy - p1,pC.fill(target_confidence),pM,1,1);
      nmask.draw_rectangle(target_x - ox - p1,target_y - oy - p1,0,0,target_x - ox + p2,target_y - oy + p2,0,0,1);
      priorities.draw_rectangle(target_x - ox - (int)patch_size,
                                target_y - oy - (int)patch_size,0,0,
                                target_x - ox + 3*p2/2,
                                target_y - oy + 3*p2/2,0,0,-1);
      // Remember patch positions.
      unsigned int *ptr_saved_patches = saved_patches.data(0,nb_saved_patches);
      *(ptr_saved_patches++) = (unsigned int)best_x;
      *(ptr_saved_patches++) = (unsigned int)best_y;
      *(ptr_saved_patches++) = (unsigned int)target_x;
      *ptr_saved_patches = (unsigned int)target_y;
      if (++nb_saved_patches>=saved_patches._height) saved_patches.resize(4,-200,1,1,0);
    }
  }
  nmask.assign(); // Free some unused memory resources
  priorities.assign();
  confidences.assign();
  is_visited.assign();

  // Blend inpainting result (if requested), using multi-scale blending algorithm.
  if (blend_size && blend_scales) {
    const float _blend_threshold = std::max(0.f,std::min(1.f,blend_threshold));
    saved_patches._height = nb_saved_patches;

    // Re-crop image and mask if outer blending is activated.
    if (is_blend_outer) {
      const int
        b2 = (int)blend_size/2, b1 = (int)blend_size - b2 - 1,
        xb0 = std::max(0,ox - b1),
        yb0 = std::max(0,oy - b1),
        xb1 = (int)std::min(_width - 1,xb0 + dx + b1 + b2),
        yb1 = (int)std::min(_height - 1,yb0 + dy + b1 + b2);
      ox = xb0; oy = yb0; dx = xb1 - xb0 + 1U, dy = yb1 - yb0 + 1U;
    }

    // Generate map of source offsets.
    CImg<unsigned int> offsets(dx,dy,1,2);
    {
      unsigned int *ptr = saved_patches.end();
      cimg_forY(saved_patches,i) {
        const unsigned int yd = *(--ptr), xd = *(--ptr), ys = *(--ptr), xs = *(--ptr);
        for (int l = -p1; l<=p2; ++l)
          for (int k = -p1; k<=p2; ++k) {
            const int xdk = (int)xd + k, ydl = (int)yd + l;
            if (xdk>=0 && xdk<=width() - 1 && ydl>=0 && ydl<=height() - 1 && mask(xd + k,yd + l)) {
              offsets((int)xd - ox + k,(int)yd - oy + l,0) = xs + k;
              offsets((int)xd - ox + k,(int)yd - oy + l,1) = ys + l;
            }
          }
      }
    }
    unsigned int *ptrx = offsets.data(0,0,0,0), *ptry = offsets.data(0,0,0,1);
    cimg_forXY(offsets,x,y) {
      if (!mask(x + ox,y + oy)) { *ptrx = (unsigned int)(x + ox); *ptry = (unsigned int)(y + oy); }
      ++ptrx; ++ptry;
    }

    // Generate map of local blending amplitudes.
    CImg<floatT> blend_map(dx,dy,1,1,0);
    CImg_3x3(I,float);
    cimg_for3XY(offsets,x,y) if (mask(x + ox,y + oy)) {
      const float
        iox = std::max((float)offsets(_n1x,y,0) - offsets(x,y,0),
                        (float)offsets(x,y,0) - offsets(_p1x,y,0)),
        ioy = std::max((float)offsets(x,_n1y,1) - offsets(x,y,1),
                        (float)offsets(x,y,1) - offsets(x,_p1y,1)),
        ion = std::sqrt(iox*iox + ioy*ioy);
      float iin = 0;
      cimg_forC(*this,c) {
        cimg_get3x3(*this,x,y,0,c,I,float);
        const float
          iix = (float)std::max(Inc - Icc,Icc - Ipc),
          iiy = (float)std::max(Icn - Icc,Icc - Icp);
        iin+=std::log(1 + iix*iix + iiy*iiy);
      }
      iin/=_spectrum;
      blend_map(x,y) = ion*iin;
    }
    blend_map.threshold(blend_map.max()*_blend_threshold).distance(1);
    cimg_forXY(blend_map,x,y) blend_map(x,y) = 1/(1 + blend_decay*blend_map(x,y));
    blend_map.quantize(blend_scales + 1,false);
    float bm, bM = blend_map.max_min(bm);
    if (bm==bM) blend_map.fill((float)blend_scales);

    // Generate blending scales.
    CImg<T> result = _inpaint_patch_crop(ox,oy,ox + dx - 1,oy + dy - 1,0);
    for (unsigned int blend_iter = 1; blend_iter<=blend_scales; ++blend_iter) {
      const unsigned int
        _blend_width = blend_iter*blend_size/blend_scales,
        blend_width = _blend_width?_blend_width + 1 - (_blend_width%2):0;
      if (!blend_width) continue;
      const int b2 = (int)blend_width/2, b1 = (int)blend_width - b2 - 1;
      CImg<floatT>
        blended = _inpaint_patch_crop(ox,oy,ox + dx - 1,oy + dy - 1,0),
        cumul(dx,dy,1,1);
      weights.assign(blend_width,blend_width,1,1,0).
        draw_gaussian((float)b1,(float)b1,blend_width/4.f,&one);
      cimg_forXY(cumul,x,y) cumul(x,y) = mask(x + ox,y + oy)?0.f:1.f;
      blended.mul(cumul);

      cimg_forY(saved_patches,l) {
        const unsigned int *ptr = saved_patches.data(0,l);
        const int
          xs = (int)*(ptr++),
          ys = (int)*(ptr++),
          xd = (int)*(ptr++),
          yd = (int)*(ptr++);
        if (xs - b1<0 || ys - b1<0 || xs + b2>=width() || ys + b2>=height()) { // Blend with partial patch
          const int
            xs0 = std::max(0,xs - b1),
            ys0 = std::max(0,ys - b1),
            xs1 = std::min(width() - 1,xs + b2),
            ys1 = std::min(height() - 1,ys + b2);
          _inpaint_patch_crop(xs0,ys0,xs1,ys1,0).move_to(pP);
          weights._inpaint_patch_crop(xs0 - xs + b1,ys0 - ys + b1,xs1 - xs + b1,ys1 - ys + b1,0).move_to(pC);
          blended.draw_image(xd + xs0 - xs - ox,yd + ys0 - ys - oy,pP,pC,-1);
          cumul.draw_image(xd + xs0 - xs - ox,yd + ys0 - ys - oy,pC,-1);
        } else { // Blend with full-size patch
          _inpaint_patch_crop(xs - b1,ys - b1,xs + b2,ys + b2,0).move_to(pP);
          blended.draw_image(xd - b1 - ox,yd - b1 - oy,pP,weights,-1);
          cumul.draw_image(xd - b1 - ox,yd - b1 - oy,weights,-1);
        }
      }

      if (is_blend_outer) {
        cimg_forXY(blended,x,y) if (blend_map(x,y)==blend_iter) {
          const float cum = cumul(x,y);
          if (cum>0) cimg_forC(*this,c) result(x,y,c) = (T)(blended(x,y,c)/cum);
        }
      } else { cimg_forXY(blended,x,y) if (mask(x + ox,y + oy) && blend_map(x,y)==blend_iter) {
          const float cum = cumul(x,y);
          if (cum>0) cimg_forC(*this,c) result(x,y,c) = (T)(blended(x,y,c)/cum);
        }
      }
    }
    if (is_blend_outer) draw_image(ox,oy,result);
    else cimg_forXY(result,x,y) if (mask(x + ox,y + oy))
           cimg_forC(*this,c) (*this)(x + ox,y + oy,c) = (T)result(x,y,c);
  }
  return *this;
}

// Special crop function that supports more boundary conditions :
// 0=dirichlet (with value 0), 1=dirichlet (with value 1) and 2=neumann.
CImg<T> _inpaint_patch_crop(const int x0, const int y0, const int x1, const int y1,
                            const unsigned int boundary=0) const {
  const int
    nx0 = x0<x1?x0:x1, nx1 = x0^x1^nx0,
    ny0 = y0<y1?y0:y1, ny1 = y0^y1^ny0;
  CImg<T> res(1U + nx1 - nx0,1U + ny1 - ny0,1,_spectrum);
  if (nx0<0 || nx1>=width() || ny0<0 || ny1>=height()) {
    if (boundary>=2) cimg_forXYZC(res,x,y,z,c) res(x,y,z,c) = _atXY(nx0 + x,ny0 + y,z,c);
    else res.fill((T)boundary).draw_image(-nx0,-ny0,*this);
  } else res.draw_image(-nx0,-ny0,*this);
  return res;
}

template<typename t>
CImg<T> get_inpaint_patch(const CImg<t>& mask, const unsigned int patch_size=11,
                          const unsigned int lookup_size=22, const float lookup_factor=1,
                          const int lookup_increment=1,
                          const unsigned int blend_size=0, const float blend_threshold=0.5,
                          const float blend_decay=0.02f, const unsigned int blend_scales=10,
                          const bool is_blend_outer=false) const {
  return (+*this).inpaint_patch(mask,patch_size,lookup_size,lookup_factor,lookup_increment,
                                blend_size,blend_threshold,blend_decay,blend_scales,is_blend_outer);
}

CImg<T>& max(const char *const expression, CImgList<T> &images) {
  return max((+*this)._fill(expression,true,3,&images,"max",this,0));
}

CImg<T>& maxabs(const char *const expression, CImgList<T> &images) {
  return maxabs((+*this)._fill(expression,true,3,&images,"maxabs",this,0));
}

CImg<T>& min(const char *const expression, CImgList<T> &images) {
  return min((+*this)._fill(expression,true,3,&images,"min",this,0));
}

CImg<T>& minabs(const char *const expression, CImgList<T> &images) {
  return minabs((+*this)._fill(expression,true,3,&images,"minabs",this,0));
}

CImg<T>& operator_andeq(const char *const expression, CImgList<T> &images) {
  return operator&=((+*this)._fill(expression,true,3,&images,"operator&=",this,0));
}

CImg<T>& operator_brseq(const char *const expression, CImgList<T> &images) {
  return operator<<=((+*this)._fill(expression,true,3,&images,"operator<<=",this,0));
}

CImg<T>& operator_blseq(const char *const expression, CImgList<T> &images) {
  return operator>>=((+*this)._fill(expression,true,3,&images,"operator>>=",this,0));
}

CImg<T>& operator_diveq(const char *const expression, CImgList<T> &images) {
  return div((+*this)._fill(expression,true,3,&images,"operator/=",this,0));
}

template<typename t>
CImg<T>& operator_eq(const t value) {
  if (is_empty()) return *this;
  cimg_openmp_for(*this,*ptr == (T)value,131072);
  return *this;
}

CImg<T>& operator_eq(const char *const expression, CImgList<T> &images) {
  return operator_eq((+*this)._fill(expression,true,3,&images,"operator_eq",this,0));
}

template<typename t>
CImg<T>& operator_eq(const CImg<t>& img) {
  const ulongT siz = size(), isiz = img.size();
  if (siz && isiz) {
    if (is_overlapped(img)) return operator_eq(+img);
    T *ptrd = _data, *const ptre = _data + siz;
    if (siz>isiz)
      for (ulongT n = siz/isiz; n; --n)
        for (const t *ptrs = img._data, *ptrs_end = ptrs + isiz; ptrs<ptrs_end; ++ptrd)
          *ptrd = (T)(*ptrd == (T)*(ptrs++));
    for (const t *ptrs = img._data; ptrd<ptre; ++ptrd) *ptrd = (T)(*ptrd == (T)*(ptrs++));
  }
  return *this;
}

template<typename t>
CImg<T>& operator_ge(const t value) {
  if (is_empty()) return *this;
  cimg_openmp_for(*this,*ptr >= (T)value,131072);
  return *this;
}

CImg<T>& operator_ge(const char *const expression, CImgList<T> &images) {
  return operator_ge((+*this)._fill(expression,true,3,&images,"operator_ge",this,0));
}

template<typename t>
CImg<T>& operator_ge(const CImg<t>& img) {
  const ulongT siz = size(), isiz = img.size();
  if (siz && isiz) {
    if (is_overlapped(img)) return operator_ge(+img);
    T *ptrd = _data, *const ptre = _data + siz;
    if (siz>isiz)
      for (ulongT n = siz/isiz; n; --n)
        for (const t *ptrs = img._data, *ptrs_end = ptrs + isiz; ptrs<ptrs_end; ++ptrd)
          *ptrd = (T)(*ptrd >= (T)*(ptrs++));
    for (const t *ptrs = img._data; ptrd<ptre; ++ptrd) *ptrd = (T)(*ptrd >= (T)*(ptrs++));
  }
  return *this;
}

template<typename t>
CImg<T>& operator_gt(const t value) {
  if (is_empty()) return *this;
  cimg_openmp_for(*this,*ptr > (T)value,131072);
  return *this;
}

CImg<T>& operator_gt(const char *const expression, CImgList<T> &images) {
  return operator_gt((+*this)._fill(expression,true,3,&images,"operator_gt",this,0));
}

template<typename t>
CImg<T>& operator_gt(const CImg<t>& img) {
  const ulongT siz = size(), isiz = img.size();
  if (siz && isiz) {
    if (is_overlapped(img)) return operator_gt(+img);
    T *ptrd = _data, *const ptre = _data + siz;
    if (siz>isiz)
      for (ulongT n = siz/isiz; n; --n)
        for (const t *ptrs = img._data, *ptrs_end = ptrs + isiz; ptrs<ptrs_end; ++ptrd)
          *ptrd = (T)(*ptrd > (T)*(ptrs++));
    for (const t *ptrs = img._data; ptrd<ptre; ++ptrd) *ptrd = (T)(*ptrd > (T)*(ptrs++));
  }
  return *this;
}

template<typename t>
CImg<T>& operator_le(const t value) {
  if (is_empty()) return *this;
  cimg_openmp_for(*this,*ptr <= (T)value,131072);
  return *this;
}

CImg<T>& operator_le(const char *const expression, CImgList<T> &images) {
  return operator_le((+*this)._fill(expression,true,3,&images,"operator_le",this,0));
}

template<typename t>
CImg<T>& operator_le(const CImg<t>& img) {
  const ulongT siz = size(), isiz = img.size();
  if (siz && isiz) {
    if (is_overlapped(img)) return operator_le(+img);
    T *ptrd = _data, *const ptre = _data + siz;
    if (siz>isiz)
      for (ulongT n = siz/isiz; n; --n)
        for (const t *ptrs = img._data, *ptrs_end = ptrs + isiz; ptrs<ptrs_end; ++ptrd)
          *ptrd = (T)(*ptrd <= (T)*(ptrs++));
    for (const t *ptrs = img._data; ptrd<ptre; ++ptrd) *ptrd = (T)(*ptrd <= (T)*(ptrs++));
  }
  return *this;
}

template<typename t>
CImg<T>& operator_lt(const t value) {
  if (is_empty()) return *this;
  cimg_openmp_for(*this,*ptr < (T)value,131072);
  return *this;
}

CImg<T>& operator_lt(const char *const expression, CImgList<T> &images) {
  return operator_lt((+*this)._fill(expression,true,3,&images,"operator_lt",this,0));
}

template<typename t>
CImg<T>& operator_lt(const CImg<t>& img) {
  const ulongT siz = size(), isiz = img.size();
  if (siz && isiz) {
    if (is_overlapped(img)) return operator_lt(+img);
    T *ptrd = _data, *const ptre = _data + siz;
    if (siz>isiz)
      for (ulongT n = siz/isiz; n; --n)
        for (const t *ptrs = img._data, *ptrs_end = ptrs + isiz; ptrs<ptrs_end; ++ptrd)
          *ptrd = (T)(*ptrd < (T)*(ptrs++));
    for (const t *ptrs = img._data; ptrd<ptre; ++ptrd) *ptrd = (T)(*ptrd < (T)*(ptrs++));
  }
  return *this;
}

CImg<T>& operator_minuseq(const char *const expression, CImgList<T> &images) {
  return operator-=((+*this)._fill(expression,true,3,&images,"operator-=",this,0));
}

CImg<T>& operator_modeq(const char *const expression, CImgList<T> &images) {
  return operator%=((+*this)._fill(expression,true,3,&images,"operator%=",this,0));
}

CImg<T>& operator_muleq(const char *const expression, CImgList<T> &images) {
  return mul((+*this)._fill(expression,true,3,&images,"operator*=",this,0));
}

template<typename t>
CImg<T>& operator_neq(const t value) {
  if (is_empty()) return *this;
  cimg_openmp_for(*this,*ptr != (T)value,131072);
  return *this;
}

CImg<T>& operator_neq(const char *const expression, CImgList<T> &images) {
  return operator_neq((+*this)._fill(expression,true,3,&images,"operator_neq",this,0));
}

template<typename t>
CImg<T>& operator_neq(const CImg<t>& img) {
  const ulongT siz = size(), isiz = img.size();
  if (siz && isiz) {
    if (is_overlapped(img)) return operator_neq(+img);
    T *ptrd = _data, *const ptre = _data + siz;
    if (siz>isiz)
      for (ulongT n = siz/isiz; n; --n)
        for (const t *ptrs = img._data, *ptrs_end = ptrs + isiz; ptrs<ptrs_end; ++ptrd)
          *ptrd = (T)(*ptrd != (T)*(ptrs++));
    for (const t *ptrs = img._data; ptrd<ptre; ++ptrd) *ptrd = (T)(*ptrd != (T)*(ptrs++));
  }
  return *this;
}

CImg<T>& operator_oreq(const char *const expression, CImgList<T> &images) {
  return operator|=((+*this)._fill(expression,true,3,&images,"operator|=",this,0));
}

CImg<T>& operator_pluseq(const char *const expression, CImgList<T> &images) {
  return operator+=((+*this)._fill(expression,true,3,&images,"operator+=",this,0));
}

CImg<T>& operator_xoreq(const char *const expression, CImgList<T> &images) {
  return operator^=((+*this)._fill(expression,true,3,&images,"operator^=",this,0));
}

CImg<T>& pow(const char *const expression, CImgList<T> &images) {
  return pow((+*this)._fill(expression,true,3,&images,"pow",this,0));
}

template<typename t>
CImg<T>& replace(CImg<t>& img) {
  return img.move_to(*this);
}

template<typename t>
CImg<T> get_replace(const CImg<t>& img) const {
  return +img;
}

CImg<T>& gmic_eval(const char *const expression, CImgList<T> &images) {
  return _fill(expression,true,6,&images,"eval",0,0);
}

CImg<T> get_gmic_eval(const char *const expression, CImgList<T> &images) const {
  return (+*this).gmic_eval(expression,images);
}

CImg<T>& gmic_fill(const char *const expression, CImgList<T> &images) {
  return _fill(expression,true,3,&images,"eval",0,0);
}

CImg<T> get_gmic_fill(const char *const expression, CImgList<T> &images) const {
  return (+*this).gmic_fill(expression,images);
}

CImg<T>& rol(const char *const expression, CImgList<T> &images) {
  return rol((+*this)._fill(expression,true,3,&images,"rol",this,0));
}

CImg<T>& ror(const char *const expression, CImgList<T> &images) {
  return ror((+*this)._fill(expression,true,3,&images,"ror",this,0));
}

template<typename t>
CImg<T>& rotate_CImg3d(const CImg<t>& rot) {
  CImg<charT> error_message(1024);
  if (!is_CImg3d(false,error_message))
    throw CImgInstanceException(_cimg_instance
                                "rotate_CImg3d(): image instance is not a CImg3d (%s).",
                                cimg_instance,error_message.data());
  const unsigned int nbv = cimg::float2uint((float)(*this)[6]);
  const T *ptrs = data() + 8;
  const float
    a = (float)rot(0,0), b = (float)rot(1,0), c = (float)rot(2,0),
    d = (float)rot(0,1), e = (float)rot(1,1), f = (float)rot(2,1),
    g = (float)rot(0,2), h = (float)rot(1,2), i = (float)rot(2,2);
  T *ptrd = data() + 8;
  for (unsigned int j = 0; j<nbv; ++j) {
    const float
      x = (float)ptrs[0],
      y = (float)ptrs[1],
      z = (float)ptrs[2];
    ptrs+=3;
    ptrd[0] = (T)(a*x + b*y + c*z);
    ptrd[1] = (T)(d*x + e*y + f*z);
    ptrd[2] = (T)(g*x + h*y + i*z);
    ptrd+=3;
  }
  return *this;
}

template<typename t>
CImg<T> get_rotate_CImg3d(const CImg<t>& rot) const {
  return (+*this).rotate_CImg3d(rot);
}

CImg<T>& scale_CImg3d(const float sx, const float sy, const float sz) {
  CImg<charT> error_message(1024);
  if (!is_CImg3d(false,error_message))
    throw CImgInstanceException(_cimg_instance
                                "scale_CImg3d(): image instance is not a CImg3d (%s).",
                                cimg_instance,error_message.data());
  const unsigned int nbv = cimg::float2uint((float)(*this)[6]);
  T *ptrd = data() + 8;
  for (unsigned int j = 0; j<nbv; ++j) { *(ptrd++)*=(T)sx; *(ptrd++)*=(T)sy; *(ptrd++)*=(T)sz; }
  return *this;
}

CImg<T> get_scale_CImg3d(const float sx, const float sy, const float sz) const {
  return (+*this).scale_CImg3d(sx,sy,sz);
}

CImg<T>& shift_CImg3d(const float tx, const float ty, const float tz) {
  CImg<charT> error_message(1024);
  if (!is_CImg3d(false,error_message))
    throw CImgInstanceException(_cimg_instance
                                "shift_CImg3d(): image instance is not a CImg3d (%s).",
                                cimg_instance,error_message.data());
  const unsigned int nbv = cimg::float2uint((float)(*this)[6]);
  T *ptrd = data() + 8;
  for (unsigned int j = 0; j<nbv; ++j) { *(ptrd++)+=(T)tx; *(ptrd++)+=(T)ty; *(ptrd++)+=(T)tz; }
  return *this;
}

CImg<T> get_shift_CImg3d(const float tx, const float ty, const float tz) const {
  return (+*this).shift_CImg3d(tx,ty,tz);
}

static const CImgList<T>& save_gmz(const char *filename, const CImgList<T>& images, const CImgList<charT>& names) {
  CImgList<T> gmz(images.size() + 1);
  cimglist_for(images,l) gmz[l].assign(images[l],true);
  CImg<charT> gmz_info = CImg<charT>::string("GMZ");
  gmz_info.append((names>'x'),'x').unroll('y').move_to(gmz.back());
  gmz.save_cimg(filename,true);
  return images;
}

//--------------- End of CImg<T> plug-in ----------------------------

// Add G'MIC-specific methods to the CImgList<T> class of the CImg library.
//-------------------------------------------------------------------------
#undef cimg_plugin
#elif defined(cimglist_plugin)

template<typename t>
static CImgList<T> copy_rounded(const CImgList<t>& list) {
  if (!cimg::type<t>::is_float() || cimg::type<T>::is_float()) return list;
  CImgList<T> res(list.size());
  cimglist_for(res,l) CImg<T>::copy_rounded(list[l]).move_to(res[l]);
  return res;
}

static CImg<T> copy_rounded(const CImg<T>& list) {
  return CImgList<T>(list,true);
}

// The method below is a variant of the method 'CImgList<T>::_display()', where
// G'MIC command 'display2d' is used in place of the built-in method 'CImg<T>::display()',
// for displaying 2d images only.
template<typename t>
const CImgList<T>& _gmic_display(CImgDisplay &disp, const char *const title, const CImgList<charT> *const titles,
                                 const bool display_info, const char axis, const float align, unsigned int *const XYZ,
                                 const bool exit_on_anykey, const unsigned int orig, const bool is_first_call,
                                 bool &is_exit,
                                 t& gmic_instance0, CImgList<T>& images, CImgList<charT>& images_names) const {
  if (is_empty())
    throw CImgInstanceException(_cimglist_instance
                                "display(): Empty instance.",
                                cimglist_instance);
  if (!disp) {
    if (axis=='x') {
      unsigned int sum_width = 0, max_height = 0;
      cimglist_for(*this,l) {
        const CImg<T> &img = _data[l];
        const unsigned int
          w = CImgDisplay::_fitscreen(img._width,img._height,img._depth,128,-85,false),
          h = CImgDisplay::_fitscreen(img._width,img._height,img._depth,128,-85,true);
        sum_width+=w;
        if (h>max_height) max_height = h;
      }
      disp.assign(cimg_fitscreen(sum_width,max_height,1U),title?title:titles?titles->__display()._data:0,1);
    } else {
      unsigned int max_width = 0, sum_height = 0;
      cimglist_for(*this,l) {
        const CImg<T> &img = _data[l];
        const unsigned int
          w = CImgDisplay::_fitscreen(img._width,img._height,img._depth,128,-85,false),
          h = CImgDisplay::_fitscreen(img._width,img._height,img._depth,128,-85,true);
        if (w>max_width) max_width = w;
        sum_height+=h;
      }
      disp.assign(cimg_fitscreen(max_width,sum_height,1U),title?title:titles?titles->__display()._data:0,1);
    }
    if (!title && !titles) disp.set_title("CImgList<%s> (%u)",pixel_type(),_width);
  } else if (title) disp.set_title("%s",title);
  else if (titles) disp.set_title("%s",titles->__display()._data);
  const CImg<char> dtitle = CImg<char>::string(disp.title());
  if (display_info) print(disp.title());
  disp.show().flush();

  if (_width==1) {
    const unsigned int dw = disp._width, dh = disp._height;
    if (!is_first_call)
      disp.resize(cimg_fitscreen(_data[0]._width,_data[0]._height,_data[0]._depth),false);
    disp.set_title("%s (%ux%ux%ux%u)",
                   dtitle.data(),_data[0]._width,_data[0]._height,_data[0]._depth,_data[0]._spectrum);
    if (_data[0]._depth==1) { // Use custom command 'display2d' for 2D images.
      CImgList<T> _images(_data[0],true);
      CImgList<charT> _images_names(dtitle,true);
      CImg<charT> com(128);
      bool is_exception = false;
      cimg_snprintf(com,com._width,"_d2d_core %d",(int)!is_first_call);
      t gmic_instance;
      cimg::swap(gmic_instance.commands,gmic_instance0.commands);
      cimg::swap(gmic_instance.commands_names,gmic_instance0.commands_names);
      cimg::swap(gmic_instance.commands_has_arguments,gmic_instance0.commands_has_arguments);
      void *const _display_window0 = gmic_instance.display_windows[0];
      gmic_instance.display_windows[0] = &disp;
      gmic_instance.is_abort = gmic_instance0.is_abort;
      try { gmic_instance.run(com.data(),_images,_images_names); }
      catch (...) { is_exception = true; }
      cimg::swap(gmic_instance.commands,gmic_instance0.commands);
      cimg::swap(gmic_instance.commands_names,gmic_instance0.commands_names);
      cimg::swap(gmic_instance.commands_has_arguments,gmic_instance0.commands_has_arguments);
      gmic_instance.display_windows[0] = _display_window0;
      if (is_exception) throw CImgDisplayException("");
    } else _data[0]._display(disp,0,false,XYZ,exit_on_anykey,!is_first_call); // Otherwise, use standard display()
    if (disp.key()) is_exit = true;
    disp.resize(cimg_fitscreen(dw,dh,1U),false).set_title("%s",dtitle.data());
  } else {
    bool disp_resize = !is_first_call;
    while (!disp.is_closed() && !is_exit) {
      const CImg<intT> s = _select(disp,0,true,axis,align,exit_on_anykey,orig,disp_resize,!is_first_call,true);
      disp_resize = true;
      if (s[0]<0 && !disp.wheel()) { // No selections done
        if (disp.button()&2) { disp.flush(); break; }
        is_exit = true;
      } else if (disp.wheel()) { // Zoom in/out
        const int wheel = disp.wheel();
        disp.set_wheel();
        if (!is_first_call && wheel<0) break;
        if (wheel>0 && _width>=4) {
          const unsigned int
            delta = std::max(1U,(unsigned int)cimg::round(0.3*_width)),
            ind0 = (unsigned int)std::max(0,s[0] - (int)delta),
            ind1 = (unsigned int)std::min(width() - 1,s[0] + (int)delta);
          if ((ind0!=0 || ind1!=_width - 1) && ind1 - ind0>=3) {
            const CImgList<T> sublist = get_shared_images(ind0,ind1);
            CImgList<charT> t_sublist;
            if (titles) t_sublist = titles->get_shared_images(ind0,ind1);
            sublist._gmic_display(disp,0,titles?&t_sublist:0,false,axis,align,XYZ,exit_on_anykey,
                                  orig + ind0,false,is_exit,
                                  gmic_instance0,images,images_names);
          }
        }
      } else if (s[0]!=0 || s[1]!=width() - 1) {
        const CImgList<T> sublist = get_shared_images(s[0],s[1]);
        CImgList<charT> t_sublist;
        if (titles) t_sublist = titles->get_shared_images(s[0],s[1]);
        sublist._gmic_display(disp,0,titles?&t_sublist:0,false,axis,align,XYZ,exit_on_anykey,
                              orig + s[0],false,is_exit,
                              gmic_instance0,images,images_names);
      }
      disp.set_title("%s",dtitle.data());
    }
  }
  return *this;
}

#undef cimglist_plugin

//--------------- End of CImgList<T> plug-in ------------------------

#else // #if defined(cimg_plugin) .. #elif defined(cimglist_plugin)

#include "gmic.h"
using namespace gmic_library;

#ifdef gmic_community
#include "gmic_stdlib_community.h"
#else
#include "gmic_stdlib.h"
#endif

// Define convenience macros, variables and functions.
//----------------------------------------------------

#undef min
#undef max

// Define number of hash slots to store variables and commands (both must be 2^n).
#ifndef gmic_varslots
#define gmic_varslots 2048
// gmic_varslots are divided in three parts:
// [ 0 -> int(gmic_varslots/2) ] : Slots for regular variables.
// [ int(gmic_varslots/2) -> int(6*gmic_varslots/7)-1 ] : Slots for global variables.
// [ int(6*gmic_varslots/7) -> gmic_varslots-1 ] : Slots for inter-thread global variables.
#endif
#ifndef gmic_comslots
#define gmic_comslots 1024
#endif
#ifndef gmic_winslots
#define gmic_winslots 10
#endif

// Macro to force stringifying selection for error messages.
#define gmic_selection_err selection2string(selection,images_names,1,gmic_selection)

inline bool is_xyzc(const char c) {
  return c=='x' || c=='y' || c=='z' || c=='c';
}

// Return image argument as a shared or non-shared copy of one existing image.
inline bool _gmic_image_arg(const unsigned int ind, const CImg<unsigned int>& selection) {
  cimg_forY(selection,l) if (selection[l]==ind) return true;
  return false;
}
#define gmic_image_arg(ind) gmic_check(_gmic_image_arg(ind,selection)?images[ind]:\
                                       images[ind].get_shared())

// Macro to manage argument substitutions from a command.
template<typename T>
void gmic::_gmic_substitute_args(const char *const argument, const char *const argument0,
                                 const char *const command, const char *const item,
                                 const CImgList<T>& images) {
  if (is_debug) {
    if (std::strcmp(argument,argument0))
      debug(images,"Command '%s': arguments = '%s' -> '%s'.",
            *command?command:item,argument0,argument);
    else
      debug(images,"Command '%s': arguments = '%s'.",
            *command?command:item,argument0);
  }
}

#define gmic_substitute_args(is_image_expr) { \
  const char *const argument0 = argument; \
  if (is_subst_arg) { \
    substitute_item(argument,images,images_names,parent_images,parent_images_names,variables_sizes,\
                    command_selection,is_image_expr).move_to(_argument); \
    argument = _argument; \
  } \
  _gmic_substitute_args(argument,argument0,command,item,images); \
}

// Macros for computing a readable version of a command argument.
inline char *_gmic_argument_text(const char *const argument, char *const argument_text, const bool is_verbose) {
  if (is_verbose) return cimg::strellipsize(argument,argument_text,80,false);
  else return &(*argument_text=0);
}
#define gmic_argument_text_printed() _gmic_argument_text(argument,gmic_use_argument_text,is_verbose)
#define gmic_argument_text() _gmic_argument_text(argument,gmic_use_argument_text,true)

// Macro for having 'get' or 'non-get' versions of G'MIC commands.
// Set 'optim_inplace' to true, only for function implementations that act 'in-place'.
#define gmic_apply(function,optim_inplace) { \
    uind = selection[l]; \
    gmic_check(images[uind]); \
    if (is_get) { \
      if (optim_inplace) \
        CImg<gmic_pixel_type>(images[uind],false).function.move_to(images); /* Surprisingly faster */ \
      else \
        images[uind].get_##function.move_to(images); \
      images_names[uind].get_copymark().move_to(images_names); \
    } else images[uind].function; \
  }

// Same as 'gmic_apply' but force computation with double-precision images.
#define gmic_apply_double(function) { \
    uind = selection[l]; \
    gmic_check(images[uind]); \
    if (is_get) { \
      CImg<double>(images[uind],false).function.move_to(images); \
      images_names[uind].get_copymark().move_to(images_names); \
    } else CImg<double>(images[uind],false).function.move_to(images[uind]); \
  }

// Macro for simple commands that has no arguments and act on images.
#define gmic_simple_command(command_name,function,description) \
  if (!std::strcmp(command_name,command)) { \
    print(images,0,description,gmic_selection.data()); \
    cimg_forY(selection,l) gmic_apply(function(),true); \
    is_change = true; \
    continue; \
}

// Macro for G'MIC arithmetic commands.
#define gmic_arithmetic_command(command_name,\
                                function1,description1,arg1_1,arg1_2,arg1_3,value_type1, \
                                function2,description2,arg2_1,arg2_2, \
                                function3,description3,arg3_1,arg3_2, \
                                description4) \
 if (!std::strcmp(command_name,command)) { \
   gmic_substitute_args(true); \
   sep = 0; value = 0; \
   if (cimg_sscanf(argument,"%lf%c",&value,&end)==1 || \
       (cimg_sscanf(argument,"%lf%c%c",&value,&sep,&end)==2 && sep=='%')) { \
     const char *const ssep = sep=='%'?"%":""; \
     print(images,0,description1 ".",arg1_1,arg1_2,arg1_3); \
     cimg_forY(selection,l) { \
       CImg<T>& img = gmic_check(images[selection[l]]); \
       nvalue = value; \
       if (sep=='%' && img) { \
         vmax = (double)img.max_min(vmin); \
         nvalue = vmin + (vmax - vmin)*value/100; \
       } \
       if (is_get) { \
         g_img.assign(img,false).function1((value_type1)nvalue).move_to(images); \
         images_names[selection[l]].get_copymark().move_to(images_names); \
       } else img.function1((value_type1)nvalue); \
     } \
     ++position; \
   } else if (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep,&end)==2 && sep==']' && \
              (ind=selection2cimg(indices,images.size(),images_names,command_name)).height()==1) { \
     print(images,0,description2 ".",arg2_1,arg2_2); \
     const CImg<T> img0 = gmic_image_arg(*ind); \
     cimg_forY(selection,l) { \
       CImg<T>& img = gmic_check(images[selection[l]]); \
       if (is_get) { \
         g_img.assign(img,false).function2(img0).move_to(images); \
         images_names[selection[l]].get_copymark().move_to(images_names); \
       } else img.function2(img0); \
     } \
     ++position; \
   } else if (cimg_sscanf(argument,"'%4095[^']%c%c",gmic_use_formula,&sep,&end)==2 && sep=='\'') { \
     strreplace_fw(formula); print(images,0,description3 ".",arg3_1,arg3_2); \
     cimg_forY(selection,l) { \
       CImg<T>& img = gmic_check(images[selection[l]]); \
       if (is_get) { \
         g_img.assign(img,false).function3((const char*)formula,images).move_to(images); \
         images_names[selection[l]].get_copymark().move_to(images_names); \
       } else img.function3((const char*)formula,images); \
     } \
     ++position; \
   } else { \
     print(images,0,description4 ".",gmic_selection.data()); \
     if (images && selection) { \
       if (is_get) { \
         g_img.assign(gmic_check(images[selection[0]]),false); \
         for (unsigned int l = 1; l<(unsigned int)selection.height(); ++l) \
           g_img.function2(gmic_check(images[selection[l]])); \
         images_names[selection[0]].get_copymark().move_to(images_names); \
         g_img.move_to(images); \
       } else if (selection.height()>=2) { \
       CImg<T>& img = gmic_check(images[selection[0]]); \
       for (unsigned int l = 1; l<(unsigned int)selection.height(); ++l) \
         img.function2(gmic_check(images[selection[l]])); \
       remove_images(images,images_names,selection,1,selection.height() - 1); \
       }}} \
   is_change = true; \
   continue; \
 }

// Return true if specified character is considered as 'blank'.
inline bool is_blank(const char x) {
  return (x>1 && x<gmic_dollar) || (x>gmic_store && x<=' ');
}

inline void _strreplace_fw(char &c) {
  switch (c) {
  case gmic_dollar : c = '$'; break;
  case gmic_lbrace : c = '{'; break;
  case gmic_rbrace : c = '}'; break;
  case gmic_comma : c = ','; break;
  case gmic_dquote : c = '\"'; break;
  }
}

inline void _strreplace_bw(char &c) {
  switch (c) {
  case '$' : c = gmic_dollar; break;
  case '{' : c = gmic_lbrace; break;
  case '}' : c = gmic_rbrace; break;
  case ',' : c = gmic_comma; break;
  case '\"' : c = gmic_dquote; break;
  }
}

// Replace special characters in a string.
char *gmic::strreplace_fw(char *const str) {
  if (str) for (char *s = str ; *s; ++s) _strreplace_fw(*s);
  return str;
}

char *gmic::strreplace_bw(char *const str) {
  if (str) for (char *s = str ; *s; ++s) _strreplace_bw(*s);
  return str;
}

//! Escape a string.
// 'res' must be a C-string large enough ('4*strlen(str) + 1' is always safe).
unsigned int gmic::strescape(const char *const str, char *const res) {
  const char *const esc = "abtnvfr";
  char *ptrd = res;
  for (const unsigned char *ptrs = (unsigned char*)str; *ptrs; ++ptrs) {
    const unsigned char c = *ptrs;
    if (c=='\\' || c=='\'' || c=='\"') { *(ptrd++) = '\\'; *(ptrd++) = c; }
    else if (c>='\a' && c<='\r') { *(ptrd++) = '\\'; *(ptrd++) = esc[c - 7]; }
    else if (c>=' ' && c<='~') *(ptrd++) = c;
    else if (c<gmic_dollar || c>gmic_dquote) {
      *(ptrd++) = '\\';
      *(ptrd++) = (char)('0' + (c>>6));
      *(ptrd++) = (char)('0' + ((c>>3)&7));
      *(ptrd++) = (char)('0' + (c&7));
    } else *(ptrd++) = c;
  }
  *ptrd = 0;
  return (unsigned int)(ptrd - res);
}

// Parse debug info string (eq. to std::sscanf(s,"%x,%x",&line_number,&file_number).
bool gmic::get_debug_info(const char *s, unsigned int &line_number, unsigned int &file_number) {
  char c = *(++s);
  bool is_digit = (c>='0' && c<='9') || (c>='a' && c<='f');
  if (is_digit) {
    unsigned int ln = 0;
    while (is_digit) {
      (ln<<=4)|=(c>='a'?c - 'a' + 10:c - '0');
      c = *(++s);
      is_digit = (c>='0' && c<='9') || (c>='a' && c<='f');
    }
    line_number = ln;
    unsigned int fn = 0;
    if (*(s++)==',') {
      c = *s;
      is_digit = (c>='0' && c<='9') || (c>='a' && c<='f');
      while (is_digit) {
        (fn<<=4)|=(c>='a'?c - 'a' + 10:c - '0');
        c = *(++s);
        is_digit = (c>='0' && c<='9') || (c>='a' && c<='f');
      }
    }
    file_number = fn;
    return true;
  }
  return false;
}

// Round a double value as with %g.
double gmic_round(const double x) {
  char tmp[32];
  double y;
  cimg_snprintf(tmp,sizeof(tmp),"%g",x);
  std::sscanf(tmp,"%lf",&y);
  return y;
}

// Manage list of all gmic runs.
inline CImgList<void*>& gmic_runs() {
  static CImgList<void*> val;
  return val;
}

inline void* get_tid() {
#if defined(__MACOSX__) || defined(__APPLE__)
  void* tid = (void*)(cimg_ulong)getpid();
#elif cimg_OS==1
#if defined(__NetBSD__) || defined(cimg_use_pthread) || cimg_display==1
  void* tid = (void*)(cimg_ulong)pthread_self();
#else
  void* tid = (void*)(cimg_ulong)syscall(SYS_gettid);
#endif
#elif cimg_OS==2
  void* tid = (void*)(cimg_ulong)GetCurrentThreadId();
#else
  void* tid = (void*)0;
#endif // #if defined(__MACOSX__) || defined(__APPLE__)
  return tid;
}

// Search G'MIC by image list (if 'p_list!=0') *or* thread_id (if 'p_list==0').
const CImg<void*> gmic::current_run(const char *const func_name, void *const p_list) {
  CImgList<void*> &grl = gmic_runs();
  void *const tid = p_list?(void*)0:get_tid();
  int p;
  for (p = grl.width() - 1; p>=0; --p) {
    const CImg<void*> &gr = grl[p];
    if (gr && ((p_list && gr[1]==(void*)p_list) || (!p_list && gr[7]==tid))) break;
  }
  if (p<0) { // Instance not found!
    if (p_list) {
      cimg::mutex(24,0);
      throw CImgArgumentException("[" cimg_appname "] Function '%s': "
                                  "Cannot determine instance of the G'MIC interpreter.",
                                  func_name);
    }
    else return CImg<void*>::empty(); // Empty instance can be returned, only when called from 'gmic_current_is_abort()'
  }
  return grl[p].get_shared(); // Return shared image
}

// Return 'is_abort' value related to current G'MIC instance.
bool* gmic::current_is_abort() {
  cimg::mutex(24);
  static bool def = false;
  CImg<void*> gr = gmic::current_run("gmic_abort_init()",0);
  bool *const res = gr?((gmic*)(gr[0]))->is_abort:&def;
  cimg::mutex(24,0);
  return res;
}

// G'MIC-related functions for the mathematical expression evaluator.
double gmic::mp_dollar(const char *const str, void *const p_list) {
  if (!(cimg::is_varname(str) ||
        ((*str=='>' || *str=='<' || *str=='!' || *str=='^' || *str=='|') && !str[1]) ||
        (*str=='{' && str[1]=='}' && !str[2])))
    throw CImgArgumentException("[" cimg_appname "_math_parser] CImg<>: Operator '$': "
                                "Invalid variable name '%s'.",
                                str);
  cimg::mutex(24);
  const CImg<void*> gr = current_run("Operator '$'",p_list);
  gmic &gmic_instance = *(gmic*)gr[0];
  CImgList<char> &images_names = *(CImgList<char>*)gr[2];
  const unsigned int *const variables_sizes = (const unsigned int*)gr[5];
  double res = cimg::type<double>::nan();

  switch (*str) {
    case '>' : case '<' :
      if (gmic_instance.nb_repeatdones || gmic_instance.nb_dowhiles || gmic_instance.nb_fordones ||
          gmic_instance.nb_foreachdones) {
        unsigned int loop_type = 0; // 0=repeatdones, 1=dowhiles, 2=fordones, 3=foreachdones
        for (int i = (int)gmic_instance.callstack.size() - 1; i>=0; --i) { // Find type of latest loop
          const CImg<char>& s = gmic_instance.callstack[i];
          if (*s=='*') {
            if (s[1]=='r') { loop_type = 0; break; }
            else if (s[1]=='d') { loop_type = 1; break; }
            else if (s[1]=='f') { loop_type = s[4]!='e'?2:3; break; }
          }
        }
        switch (loop_type) {
        case 0 : { // repeat...done
          const unsigned int *const rd = gmic_instance.repeatdones.data(0,gmic_instance.nb_repeatdones - 1);
          res = (double)(*str=='>'?rd[1]:rd[2] - 1);
        } break;
        case 1 : { // do...while
          const unsigned int *const dw = gmic_instance.dowhiles.data(0,gmic_instance.nb_dowhiles - 1);
          if (*str=='>') res = (double)dw[1];
        } break;
        case 2 : { // for...done
          const unsigned int *const fd = gmic_instance.fordones.data(0,gmic_instance.nb_fordones - 1);
          if (*str=='>') res = (double)fd[1];
        } break;
        case 3 : { // foreach...done
          const unsigned int *const fed = gmic_instance.foreachdones.data(0,gmic_instance.nb_foreachdones - 1);
          res = (double)(*str=='>'?fed[0]:fed[1] - 1);
        } break;
        }
      }
      break;
  case '!' :
    res = (double)images_names.size();
    break;
  case '^' :
    res = (double)gmic_instance.verbosity;
    break;
  case '|' :
    res = (cimg::time() - gmic_instance.reference_time)/1000.;
    break;
  default : {
    const CImg<char> value = *str=='{'?gmic_instance.status.get_shared():
      gmic_instance.get_variable(str,variables_sizes,&images_names);
    if (value && *value) {
      char end;
      if (std::sscanf(value,"%lf%c",&res,&end)!=1) res = 0;
    }
  }
  }
  cimg::mutex(24,0);
  return res;
}

double gmic::mp_abort() {
#if defined(cimg_use_abort) && !defined(__MACOSX__) && !defined(__APPLE__)
  cimg_abort_init;
  *gmic_is_abort = true;
  cimg_abort_test;
#endif
  return cimg::type<double>::nan();
}

template<typename T>
double gmic::mp_get(double *const ptrd, const unsigned int siz, const bool to_string, const char *const str,
                    void *const p_list, const T& pixel_type) {
  cimg::unused(pixel_type);
  cimg::mutex(24);
  const CImg<void*> gr = current_run("Function 'get()'",p_list);
  gmic &gmic_instance = *(gmic*)gr[0];
  CImgList<char>& images_names = *(CImgList<char>*)gr[2];
  const unsigned int *const variables_sizes = (const unsigned int*)gr[5];
  CImg<char> _varname(256);
  char *const varname = _varname.data(), end;

  if ((cimg_sscanf(str,"%255[a-zA-Z0-9_]%c",&(*varname=0),&end)==1 && (*varname<'0' || *varname>'9')) ||
      (*str=='{' && str[1]=='}' && !str[2])) {
    const CImg<char> value = *str=='{'?gmic_instance.status.get_shared():
      gmic_instance.get_variable(varname,variables_sizes,&images_names);

    if (!value) { // Undefined variable
      if (!siz) *ptrd = cimg::type<double>::nan();
      else for (unsigned int i = 0; i<siz; ++i) ptrd[i] = cimg::type<double>::nan();
    } else if (to_string) { // Return variable content as a string
      if (!siz) { char c = *value; _strreplace_fw(c); *ptrd = c; }
      else {
        CImg<double> dest(ptrd,siz,1,1,1,true);
        CImg<char> _value(value,false);
        strreplace_fw(_value);
        dest.draw_image(_value);
        if (dest.width()>_value.width()) dest.get_shared_points(_value.width(),dest.width() - 1).fill(0);
      }
    } else { // Convert variable content as numbers
      double dvalue = 0;
      if (!siz) { // Scalar result
        if (cimg_sscanf(value,"%lf",&dvalue)!=1) *ptrd = cimg::type<double>::nan();
        else *ptrd = dvalue;
      } else { // Vector result
        CImg<double> dest(ptrd,siz,1,1,1,true);
        if (*value==gmic_store) { // Image-encoded variable
          const char *const zero = (char*)::std::memchr(value,0,value.size());
          CImgList<T> list;
          if (zero) CImgList<T>::get_unserialize(value,zero + 1 - value.data()).move_to(list);
          if (list.size()!=2) {
            cimg::mutex(24,0);
            throw CImgArgumentException("[" cimg_appname "_math_parser] CImg<%s>: Function 'get()': "
                                        "Variable '%s' stores %u images, cannot be returned as a single vector.",
                                        cimg::type<T>::string(),str,list.size());
          }
          dest = list[0].resize(siz,1,1,1,-1);

        } else { // Regular string variable
          if (cimg_sscanf(value,"%lf%c",&dvalue,&end)==1) {
            dest[0] = dvalue;
            if (dest._width>1) dest.get_shared_points(1,dest._width - 1).fill(0);
          } else if (dest.fill(0)._fill_from_values(value,false))
            for (unsigned int i = 0; i<siz; ++i) ptrd[i] = cimg::type<double>::nan();
        }
      }
    }
  } else {
    cimg::mutex(24,0);
    throw CImgArgumentException("[" cimg_appname "_math_parser] CImg<%s>: Function 'get()': "
                                "Invalid variable name '%s'.",
                                cimg::type<T>::string(),str);
  }
  cimg::mutex(24,0);
  return siz?cimg::type<double>::nan():*ptrd;
}

double gmic::mp_set(const double *const ptrs, const unsigned int siz, const char *const str,
                    void *const p_list) {
  cimg::mutex(24);
  const CImg<void*> gr = current_run("Function 'set()'",p_list);
  gmic &gmic_instance = *(gmic*)gr[0];
  const unsigned int *const variables_sizes = (const unsigned int*)gr[5];
  CImg<char> _varname(256);
  char *const varname = _varname.data(), end;

  if ((cimg_sscanf(str,"%255[a-zA-Z0-9_]%c",&(*varname=0),&end)==1 && (*varname<'0' || *varname>'9')) ||
      (*str=='{' && str[1]=='}' && !str[2])) {
    CImg<char> s_value;
    if (siz) { // Value is a string
      s_value.assign(siz + 1);
      cimg_for_inX(s_value,0,s_value.width() - 1,i) s_value[i] = (char)ptrs[i];
      s_value.back() = 0;
    } else { // Value is a scalar
      s_value.assign(24);
      cimg_snprintf(s_value,s_value.width(),"%.17g",*ptrs);
    }
    if (*str=='{') CImg<char>::string(s_value).move_to(gmic_instance.status);
    else gmic_instance.set_variable(str,'=',s_value,0,variables_sizes);
  } else {
    cimg::mutex(24,0);
    throw CImgArgumentException("[" cimg_appname "_math_parser] CImg<>: Function 'set()': "
                                "Invalid variable name '%s'.",
                                str);
  }
  cimg::mutex(24,0);
  return siz?cimg::type<double>::nan():*ptrs;
}

double gmic::mp_name(const unsigned int ind, double *const out_str, const unsigned int siz,
                     void *const p_list) {
  cimg::mutex(24);
  const CImg<void*> gr = current_run("Function 'name()'",p_list);
  CImgList<char> &images_names = *(CImgList<char>*)gr[2];
  std::memset(out_str,0,siz*sizeof(double));
  if (ind<images_names.size()) {
    const char *ptrs = images_names[ind];
    unsigned int k;
    for (k = 0; k<siz && ptrs[k]; ++k) out_str[k] = (double)ptrs[k];
    if (k<siz) out_str[k] = 0;
  }
  cimg::mutex(24,0);
  return cimg::type<double>::nan();
}

// This method is not thread-safe. Ensure it's never run in parallel!
template<typename T>
double gmic::mp_run(char *const str,
                    void *const p_list, const T& pixel_type) {
  cimg::unused(pixel_type);
  const CImg<void*> gr = current_run("Function 'run()'",p_list);
  double res = cimg::type<double>::nan();

  gmic &gmic_instance = *(gmic*)gr[0];
  CImgList<T> &images = *(CImgList<T>*)gr[1];
  CImgList<char> &images_names = *(CImgList<char>*)gr[2];
  CImgList<T> &parent_images = *(CImgList<T>*)gr[3];
  CImgList<char> &parent_images_names = *(CImgList<char>*)gr[4];
  const unsigned int *const variables_sizes = (const unsigned int*)gr[5];
  const CImg<unsigned int> *const command_selection = (const CImg<unsigned int>*)gr[6];

  CImg<char> is_error;
  char sep;
  if (gmic_instance.is_debug_info && gmic_instance.debug_line!=~0U) {
    CImg<char> title(32);
    cimg_snprintf(title,title.width(),"*expr#%u",gmic_instance.debug_line);
    CImg<char>::string(title).move_to(gmic_instance.callstack);
  } else CImg<char>::string("*expr").move_to(gmic_instance.callstack);
  unsigned int pos = 0;
  try {
    gmic_instance._run(gmic_instance.commands_line_to_CImgList(gmic::strreplace_fw(str)),pos,images,images_names,
                       parent_images,parent_images_names,variables_sizes,0,0,command_selection,false);
  } catch (gmic_exception &e) {
    CImg<char>::string(e.what()).move_to(is_error);
  }
  gmic_instance.callstack.remove();
  if (is_error || !gmic_instance.status || !*gmic_instance.status ||
      cimg_sscanf(gmic_instance.status,"%lf%c",&res,&sep)!=1)
    res = cimg::type<double>::nan();
  if (is_error)
    throw CImgArgumentException("[" cimg_appname "_math_parser] CImg<%s>: Function 'run()': %s",
                                cimg::type<T>::string(),is_error.data());
  return res;
}

template<typename T>
double gmic::mp_store(const double *const ptrs, const unsigned int siz,
                      const unsigned int w, const unsigned int h, const unsigned d, const unsigned int s,
                      const bool is_compressed, const char *const str,
                      void *const p_list, const T& pixel_type) {
  cimg::unused(pixel_type);
  cimg::mutex(24);
  const CImg<void*> gr = current_run("Function 'store()'",p_list);
  cimg_pragma_openmp(critical(mp_store))
  {
    gmic &gmic_instance = *(gmic*)gr[0];
    const unsigned int *const variables_sizes = (const unsigned int*)gr[5];
    CImg<char> _varname(256);
    char *const varname = _varname.data(), end;

    if (cimg_sscanf(str,"%255[a-zA-Z0-9_]%c",&(*varname=0),&end)==1 &&
        (*varname<'0' || *varname>'9')) {
      CImgList<T> g_list;
      const unsigned int rsiz = w*h*d*s;
      if (rsiz<=siz) CImg<T>(ptrs,w,h,d,s).move_to(g_list);
      else CImg<T>(ptrs,siz,1,1,1).resize(w,h,d,s,-1).move_to(g_list);

      CImg<char> name = CImg<char>::string(varname);
      name.resize(name.width() + 4,1,1,1,0,0,1);
      name[0] = 'G'; name[1] = 'M'; name[2] = 'Z'; name[3] = 0;
      name.unroll('y').move_to(g_list);
      g_list.get_serialize(is_compressed,(unsigned int)(9 + std::strlen(varname))).move_to(name);
      cimg_snprintf(name,name._height,"%c*store/%s",gmic_store,_varname.data());
      gmic_instance.set_variable(_varname.data(),name,variables_sizes);
    } else {
      cimg::mutex(24,0);
      throw CImgArgumentException("[" cimg_appname "_math_parser] CImg<%s>: Function 'store()': "
                                  "Invalid variable name '%s'.",
                                  cimg::type<T>::string(),str);
    }
  }
  cimg::mutex(24,0);
  return cimg::type<double>::nan();
}

// Manage mutexes.
struct _gmic_mutex {
#if cimg_OS==1 && (defined(cimg_use_pthread) || cimg_display==1)
  pthread_mutex_t mutex[256];
  _gmic_mutex() { for (unsigned int i = 0; i<256; ++i) pthread_mutex_init(&mutex[i],0); }
  void lock(const unsigned int n) { pthread_mutex_lock(&mutex[n]); }
  void unlock(const unsigned int n) { pthread_mutex_unlock(&mutex[n]); }
#elif cimg_OS==2 // #if cimg_OS==1 && (defined(cimg_use_pthread) || cimg_display==1)
  HANDLE mutex[256];
  _gmic_mutex() { for (unsigned int i = 0; i<256; ++i) mutex[i] = CreateMutex(0,FALSE,0); }
  void lock(const unsigned int n) { WaitForSingleObject(mutex[n],INFINITE); }
  void unlock(const unsigned int n) { ReleaseMutex(mutex[n]); }
#else // #if cimg_OS==1 && (defined(cimg_use_pthread) || cimg_display==1)
  _gmic_mutex() {}
  void lock(const unsigned int) {}
  void unlock(const unsigned int) {}
#endif // #if cimg_OS==1 && (defined(cimg_use_pthread) || cimg_display==1)
};
inline _gmic_mutex& gmic_mutex() { static _gmic_mutex val; return val; }

// Thread structure and routine for command 'parallel'.
template<typename T>
struct _gmic_parallel {
  CImgList<char> *images_names, *parent_images_names, commands_line;
  CImg<_gmic_parallel<T> > *gmic_threads;
  CImgList<T> *images, *parent_images;
  CImg<unsigned int> variables_sizes;
  const CImg<unsigned int> *command_selection;
  bool is_thread_running;
  gmic_exception exception;
  gmic gmic_instance;
#ifdef gmic_is_parallel
#ifdef PTHREAD_CANCEL_ENABLE
  pthread_t thread_id;
#elif cimg_OS==2
  HANDLE thread_id;
#endif // #ifdef PTHREAD_CANCEL_ENABLE
#endif // #ifdef gmic_is_parallel
  _gmic_parallel() { variables_sizes.assign(gmic_varslots); }
};

template<typename T>
#if cimg_OS==2 && !defined(PTHREAD_CANCEL_ENABLE)
DWORD WINAPI gmic_parallel(LPVOID arg) {
#else
static void *gmic_parallel(void *arg) {
#endif
  _gmic_parallel<T> &st = *(_gmic_parallel<T>*)arg;
  try {
    unsigned int pos = 0;
    st.gmic_instance.is_debug_info = false;
    st.gmic_instance._run(st.commands_line,pos,*st.images,*st.images_names,
                          *st.parent_images,*st.parent_images_names,
                          st.variables_sizes,0,0,st.command_selection,true);
  } catch (gmic_exception &e) {
    cimg_forY(*st.gmic_threads,l)
      (*st.gmic_threads)[l].gmic_instance.is_abort_thread = true;
    st.exception._command.assign(e._command);
    st.exception._message.assign(e._message);
  }
#if defined(gmic_is_parallel) && defined(PTHREAD_CANCEL_ENABLE)
  pthread_exit(0);
#endif // #if defined(gmic_is_parallel) && defined(PTHREAD_CANCEL_ENABLE)
  return 0;
}

// Array of G'MIC built-in commands (must be sorted in lexicographic order!).
const char *gmic::builtin_commands_names[] = {

  // Commands of length>3.
  "acos","acosh","add3d","append","asin","asinh","atan","atan2","atanh","autocrop",
  "bilateral","blur","boxfilter","break",
  "camera","check","check3d","command","continue","convolve","correlate","cosh","crop","cumulate","cursor",
  "debug","delete","denoise","deriche","dijkstra","dilate","discard","displacement","display","distance","div3d","done",
  "echo","eigen","eikonal","elif","ellipse","else","endian","equalize","erode","error","eval","exec",
  "files","fill","flood","foreach",
  "graph","guided",
  "histogram",
  "ifft","image","index","inpaint","input","invert","isoline3d","isosurface3d",
  "keep",
  "label","light3d","line","local","log10","log2",
  "mandelbrot","matchpatch","maxabs","mdiv","median","minabs","mirror","mmul","move","mproj","mul3d","mutex",
  "name","named","network","noarg","noise","normalize",
  "object3d","onfail","output",
  "parallel","pass","permute","plasma","plot","point","polygon","print","progress",
  "quit",
  "rand","remove","repeat","resize","return","reverse","rotate","rotate3d","round",
  "screen","select","serialize","shared","shift","sign","sinc","sinh","skip",
    "smooth","solve","sort","split","sqrt","srand","status","store","streamline3d","sub3d",
  "tanh","text","trisolve",
  "uncommand","unroll","unserialize",
  "vanvliet","verbose",
  "wait","warn","warp","watershed","while","window",
  0,

  // Commands of length 3.
  "*3d","+3d","-3d","/3d","abs","add","and","bsl","bsr","cos","cut","div","erf","exp","fft","for","j3d","l3d","log",
  "map","max","min","mod","mul","neq","nmd","pow","r3d","rol","ror","set","sin","sqr","sub","svd","tan","xor",

  // Commands of length 2.
  "!=","<<","<=","==","=>",">=",">>",
  "do","eq","fi","ge","gt","if","le","lt","m*","m/","mv","nm","or","rm","rv","sh","um",
  "w0","w1","w2","w3","w4","w5","w6","w7","w8","w9",

  // Commands of length 1.
  "%","&","*","+","-","/","<","=",">",
  "a","b","c","d","e","f","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z",
  "y","z","^","{","|","}"
};

CImg<int> gmic::builtin_commands_inds = CImg<int>::empty();

// Perform a dichotomic search in a lexicographic ordered 'CImgList<char>' or 'char**'.
// Return false or true if search succeeded.
template<typename T>
bool gmic::search_sorted(const char *const str, const T& list, const unsigned int length, unsigned int &out_ind) {
  if (!length) { out_ind = 0; return false; }
  int err, pos, posm = 0, posM = length - 1;
  do {
    pos = (posm + posM)/2;
    err = std::strcmp(list[pos],str);
    if (!err) { posm = pos; break; }
    if (err>0) posM = pos - 1; else posm = pos + 1;
  } while (posm<=posM);
  out_ind = posm;
  return !err;
}

// Return Levenshtein distance between two strings.
// (adapted from http://rosettacode.org/wiki/Levenshtein_distance#C)
int gmic::_levenshtein(const char *const s, const char *const t,
                       CImg<int>& d, const int i, const int j) {
  const int ls = d.width() - 1, lt = d.height() - 1;
  if (d(i,j)>=0) return d(i,j);
  int x;
  if (i==ls) x = lt - j;
  else if (j==lt) x = ls - i;
  else if (s[i]==t[j]) x = _levenshtein(s,t,d,i + 1,j + 1);
  else {
    x = _levenshtein(s,t,d,i + 1,j + 1);
    int y;
    if ((y=_levenshtein(s,t,d,i,j + 1))<x) x = y;
    if ((y=_levenshtein(s,t,d,i + 1,j))<x) x = y;
    ++x;
  }
  return d(i,j) = x;
}

int gmic::levenshtein(const char *const s, const char *const t) {
  const char *const ns = s?s:"", *const nt = t?t:"";
  const int ls = (int)std::strlen(ns), lt = (int)std::strlen(nt);
  if (!ls) return lt; else if (!lt) return ls;
  CImg<int> d(1 + ls,1 + lt,1,1,-1);
  return _levenshtein(ns,nt,d,0,0);
}

// Wait for threads to finish.
template<typename T>
void gmic::wait_threads(void *const p_gmic_threads, const bool try_abort, const T& pixel_type) {
  cimg::unused(p_gmic_threads,try_abort,pixel_type);
#ifdef gmic_is_parallel
  CImg<_gmic_parallel<T> > &gmic_threads = *(CImg<_gmic_parallel<T> >*)p_gmic_threads;
  cimg_forY(gmic_threads,l) {
    if (try_abort && gmic_threads[l].is_thread_running)
      gmic_threads[l].gmic_instance.is_abort_thread = true;

    cimg::mutex(25);
    if (gmic_threads[l].is_thread_running) {
      gmic_threads[l].is_thread_running = false;
      cimg::mutex(25,0);
#ifdef PTHREAD_CANCEL_ENABLE
      pthread_join(gmic_threads[l].thread_id,0);
#elif cimg_OS==2 // #ifdef PTHREAD_CANCEL_ENABLE
      WaitForSingleObject(gmic_threads[l].thread_id,INFINITE);
      CloseHandle(gmic_threads[l].thread_id);
#endif // #ifdef PTHREAD_CANCEL_ENABLE
    } else cimg::mutex(25,0);

    is_change|=gmic_threads[l].gmic_instance.is_change;
  }
#endif // #ifdef gmic_is_parallel
}

// Return a hashcode from a string.
unsigned int gmic::hashcode(const char *const str, const bool is_variable) {
  if (!str) return 0U;
  unsigned int hash = 5381U;
  if (is_variable) {
    int k = 0; for (const char *s = str; *s && k++<32; ++s) (hash*=31)+=*s;
    if (*str=='_') {
      if (str[1]=='_') return 6*gmic_varslots/7 + (hash%(gmic_varslots - 6*gmic_varslots/7));
      return gmic_varslots/2 + (hash%(6*gmic_varslots/7 - gmic_varslots/2));
    }
    return hash&(gmic_varslots/2 - 1);
  }
  int k = 0; for (const char *s = str; *s && k++<32; ++s) (hash*=31)+=*s;
  return hash&(gmic_comslots - 1);
}

// Tells if the implementation of a G'MIC command contains arguments.
bool gmic::command_has_arguments(const char *const command) {
  if (!command || !*command) return false;
  for (const char *s = std::strchr(command,'$'); s; s = std::strchr(s,'$')) {
    const char c = *(++s);
    if (c=='#' ||
        c=='*' ||
        c=='=' ||
        (c>'0' && c<='9') ||
        (c=='-' && *(s + 1)>'0' && *(s + 1)<='9') ||
        (c=='\"' && *(s + 1)=='*' && *(s + 2)=='\"') ||
        (c=='{' && (*(s + 1)=='^' ||
                    (*(s + 1)>'0' && *(s + 1)<='9') ||
                    (*(s + 1)=='-' && *(s + 2)>'0' && *(s + 2)<='9')))) return true;
  }
  return false;
}

// Compute the basename of a filename.
const char* gmic::basename(const char *const str)  {
  if (!str || !*str) return "";
  const unsigned int l = (unsigned int)std::strlen(str);
  unsigned int ll = l - 1; // 'Last' character to check
  while (ll>=3 && str[ll]>='0' && str[ll]<='9') --ll;
  if (ll>=3 && ll!=l - 1 && str[ll - 1]=='_' && str[ll]=='c' && str[ll + 1]!='0') ll-=2; // Ignore copy mark
  else ll = l - 1;
  if (*str=='[' && (str[ll]==']' || str[ll]=='.')) return str;
  const char *p = 0, *np = str;
  while (np>=str && (p=np)) np = std::strchr(np,'/') + 1;
  np = p;
  while (np>=str && (p=np)) np = std::strchr(np,'\\') + 1;
  return p;
}

// Constructors / destructors.
//----------------------------
#define gmic_display_window(n) (*(CImgDisplay*)display_windows[n])

CImg<char> gmic::stdlib = CImg<char>::empty();

gmic::gmic():gmic_new_attr {
  assign();
}

gmic& gmic::assign() {
  CImgList<gmic_pixel_type> images;
  CImgList<char> images_names;
  return _gmic(0,images,images_names,0,true,0,0);
}

template<typename T>
gmic::gmic(const char *const commands_line, const char *const custom_commands,
           const bool include_stdlib, float *const p_progress, bool *const p_is_abort,
           const T& pixel_type):gmic_new_attr {
  assign(commands_line,custom_commands,include_stdlib,p_progress,p_is_abort,pixel_type);
}

template<typename T>
gmic::gmic(const char *const commands_line, CImgList<T>& images, CImgList<char>& images_names,
           const char *const custom_commands, const bool include_stdlib,
           float *const p_progress, bool *const p_is_abort):gmic_new_attr {
  assign(commands_line,images,images_names,custom_commands,include_stdlib,p_progress,p_is_abort);
}

template<typename T>
gmic& gmic::assign(const char *const commands_line, const char *const custom_commands,
                   const bool include_stdlib, float *const p_progress, bool *const p_is_abort,
                   const T& pixel_type) {
  cimg::unused(pixel_type);
  CImgList<T> images;
  CImgList<char> images_names;
  return _gmic(commands_line,
               images,images_names,custom_commands,
               include_stdlib,p_progress,p_is_abort);
}

template<typename T>
gmic& gmic::assign(const char *const commands_line, CImgList<T>& images, CImgList<char>& images_names,
             const char *const custom_commands, const bool include_stdlib,
             float *const p_progress, bool *const p_is_abort) {
  return _gmic(commands_line,images,images_names,custom_commands,include_stdlib,p_progress,p_is_abort);
}

gmic::~gmic() {
  cimg_forX(display_windows,l) delete &gmic_display_window(l);
  delete[] commands;
  delete[] commands_names;
  delete[] commands_has_arguments;
  delete[] _variables;
  delete[] _variables_names;
  delete[] _variables_lengths;
  delete[] variables;
  delete[] variables_names;
  delete[] variables_lengths;
  cimg::exception_mode(cimg_exception_mode);
}

// Decompress G'MIC standard library commands.
//---------------------------------------------
const CImg<char>& gmic::decompress_stdlib() {
  cimg::mutex(22);
  if (!stdlib) try {
      CImgList<char>::get_unserialize(CImg<unsigned char>(data_gmic,1,size_data_gmic,1,1,true))[0].
        move_to(stdlib);
    } catch (...) {
      cimg::mutex(29);
      std::fprintf(cimg::output(),
                   "[gmic] %s*** Warning *** Could not decompress G'MIC standard library, ignoring it.%s\n",
                   cimg::t_red,cimg::t_normal);
      std::fflush(cimg::output());
      cimg::mutex(29,0);
      stdlib.assign(1,1,1,1,0);
    }
  cimg::mutex(22,0);
  return stdlib;
}

// Get the value of an environment variable.
//-------------------------------------------
static const char* gmic_getenv(const char *const varname) {
#if cimg_OS==2
  static CImg<char> utf8Buffer(768);
  // Get the value of the environment variable using the wide-character
  // (UTF-16) version of the Windows API and convert it to UTF-8.

  // Convert the environment variable name to a wide string.
  const int wideNameLength = MultiByteToWideChar(CP_UTF8,0,varname,-1,0,0);
  if (wideNameLength) {
    CImg<wchar_t> wideVarName(wideNameLength);
    if (MultiByteToWideChar(CP_UTF8,0,varname,-1,wideVarName,wideNameLength)) {
      const DWORD wideValueLength = GetEnvironmentVariableW(wideVarName,0,0);
      if (wideValueLength) {
        CImg<wchar_t> wideValue(wideValueLength);
        if (GetEnvironmentVariableW(wideVarName,wideValue,wideValueLength)) {
          // Convert the returned value from UTF-16 to UTF-8.
          const int utf8Length = WideCharToMultiByte(CP_UTF8,0,wideValue,wideValueLength,0,0,0,0);
          if (utf8Length && utf8Length<utf8Buffer.width()) {
            if (WideCharToMultiByte(CP_UTF8,0,wideValue,wideValueLength,utf8Buffer,utf8Length,0,0)) return utf8Buffer;
          }
        }
      }
    }
  }
#endif // cimgOS==2
  return getenv(varname);
}

// Get path to .gmic user file.
//-----------------------------
const char* gmic::path_user(const char *const custom_path) {
  static CImg<char> path_user;
  if (path_user) return path_user;
  cimg::mutex(28);
  const char *_path_user = 0;
  if (custom_path && cimg::is_directory(custom_path)) _path_user = custom_path;
  if (!_path_user) _path_user = gmic_getenv("GMIC_PATH");
#if cimg_OS!=2
  if (!_path_user) _path_user = gmic_getenv("HOME");
#else
  if (!_path_user) _path_user = gmic_getenv("USERPROFILE");
#endif
  if (!_path_user) _path_user = gmic_getenv("TMP");
  if (!_path_user) _path_user = gmic_getenv("TEMP");
  if (!_path_user) _path_user = gmic_getenv("TMPDIR");
  if (!_path_user) _path_user = "";
  path_user.assign(1024);
#if cimg_OS!=2
  cimg_snprintf(path_user,path_user.width(),"%s%c.gmic",
                _path_user,cimg_file_separator);
#else
  cimg_snprintf(path_user,path_user.width(),"%s%cuser.gmic",
                _path_user,cimg_file_separator);
#endif
  CImg<char>::string(path_user).move_to(path_user); // Optimize length
  cimg::mutex(28,0);
  return path_user;
}

// Get path to the resource directory.
//------------------------------------
const char* gmic::path_rc(const char *const custom_path) {
  static CImg<char> path_rc;
  CImg<char> path_tmp;
  if (path_rc) return path_rc;
  cimg::mutex(28);
  const char *_path_rc = 0;
  bool add_gmic_subfolder = true;
  if (custom_path && cimg::is_directory(custom_path)) { _path_rc = custom_path; add_gmic_subfolder = false; }
  if (!_path_rc) { _path_rc = gmic_getenv("GMIC_PATH"); add_gmic_subfolder = _path_rc==0; }
  if (!_path_rc) _path_rc = gmic_getenv("XDG_CONFIG_HOME");
#if cimg_OS!=2
  if (!_path_rc) {
    _path_rc = gmic_getenv("HOME");
    if (_path_rc) {
      path_tmp.assign(std::strlen(_path_rc) + 10);
      cimg_snprintf(path_tmp,path_tmp._width,"%s/.config",_path_rc);
      if (cimg::is_directory(path_tmp)) _path_rc = path_tmp;
    }
  }
#else
  if (!_path_rc) _path_rc = gmic_getenv("APPDATA");
#endif
  if (!_path_rc) _path_rc = gmic_getenv("TMP");
  if (!_path_rc) _path_rc = gmic_getenv("TEMP");
  if (!_path_rc) _path_rc = gmic_getenv("TMPDIR");
  if (!_path_rc) _path_rc = "";
  path_rc.assign(1024);

  if (add_gmic_subfolder)
    cimg_snprintf(path_rc,path_rc.width(),"%s%cgmic%c",_path_rc,cimg_file_separator,cimg_file_separator);
  else
    cimg_snprintf(path_rc,path_rc.width(),"%s%c",_path_rc,cimg_file_separator);

  CImg<char>::string(path_rc).move_to(path_rc); // Optimize length
  cimg::mutex(28,0);
  return path_rc;
}

// Create resources directory.
//----------------------------
bool gmic::init_rc(const char *const custom_path) {
  CImg<char> dirname = CImg<char>::string(path_rc(custom_path));
  if (dirname.width()>=2) {
    char &c = dirname[dirname.width() - 2];
    if (c=='/' || c=='\\') c = 0;
  }
  if (!cimg::is_directory(dirname)) {
#if cimg_OS==2
    DeleteFileA(dirname); // In case 'dirname' is already a file
    if (!CreateDirectoryA(dirname,0)) {
      // The path may be UTF-8, convert it to a
      // wide-character string and try again.
      const int wideLength = MultiByteToWideChar(CP_UTF8,0,dirname,-1,0,0);
      if (!wideLength) return false;
      CImg<wchar_t> wpath(wideLength);
      if (!MultiByteToWideChar(CP_UTF8,0,dirname,-1,wpath,wideLength)) return false;
      DeleteFileW(wpath);
      return (bool)CreateDirectoryW(wpath,0);
    }
#else
    std::remove(dirname); // In case 'dirname' is already a file
    return !(bool)mkdir(dirname,0777);
#endif
  }
  return true;
}

// Get current call stack as a string.
//-------------------------------------
CImg<char> gmic::callstack2string(const CImg<unsigned int> *const callstack_selection, const bool _is_debug) const {
  if (callstack_selection && !*callstack_selection) return CImg<char>("./",3);
  CImgList<char> input_callstack;
  if (!callstack_selection) input_callstack.assign(callstack,true);
  else cimg_forY(*callstack_selection,l) input_callstack.insert(callstack[(*callstack_selection)[l]],~0U,true);
  CImgList<char> res;
  const unsigned int siz = (unsigned int)input_callstack.size();
  if (siz<=9 || _is_debug) res.assign(input_callstack,false);
  else {
    res.assign(9);
    res[0].assign(input_callstack[0],false);
    res[1].assign(input_callstack[1],false);
    res[2].assign(input_callstack[2],false);
    res[3].assign(input_callstack[3],false);
    res[4].assign("(...)",6);
    res[5].assign(input_callstack[siz - 4],false);
    res[6].assign(input_callstack[siz - 3],false);
    res[7].assign(input_callstack[siz - 2],false);
    res[8].assign(input_callstack[siz - 1],false);
  }
  cimglist_for(res,l) if (res(l,0)) res[l].back() = '/'; else res.remove(l--);
  CImg<char>::vector(0).move_to(res);
  return res>'x';
}

CImg<char> gmic::callstack2string(const bool _is_debug) const {
  return callstack2string(0,_is_debug);
}

CImg<char> gmic::callstack2string(const CImg<unsigned int>& callstack_selection, const bool _is_debug) const {
  return callstack2string(&callstack_selection,_is_debug);
}

// Pop callstack until it reaches a certain size.
//-----------------------------------------------
// Used to ensure that callstack stays coherent when errors occurs in '_run()'.
void gmic::pop_callstack(const unsigned int callstack_size) {
  while (callstack.size()>callstack_size) {
    const char *const s = callstack.back();
    if (*s=='*') switch (s[1]) {
      case 'r' : --nb_repeatdones; break;
      case 'd' : --nb_dowhiles; break;
      case 'f' : if (s[4]!='e') --nb_fordones; else --nb_foreachdones; break;
      }
    callstack.remove();
  }
}

// Parse items from a G'MIC command line.
//---------------------------------------
CImgList<char> gmic::commands_line_to_CImgList(const char *const commands_line) {
  if (!commands_line || !*commands_line) return CImgList<char>();
  bool is_dquoted = false, is_subst = false;
  const char *ptrs0 = commands_line;
  while (is_blank(*ptrs0)) ++ptrs0; // Remove leading spaces to first item
  CImg<char> item((unsigned int)std::strlen(ptrs0) + 2);
  CImgList<char> items;
  char *ptrd = item.data(), c = 0;

  for (const char *ptrs = ptrs0; *ptrs; ++ptrs) {
    c = *ptrs;
    if (c=='\\') { // If escaped character
      c = *(++ptrs);
      switch (c) {
      case 0 : c = '\\'; --ptrs; break;
      case '$' : c = gmic_dollar; break;
      case '{' : c = gmic_lbrace; break;
      case '}' : c = gmic_rbrace; break;
      case ',' : c = gmic_comma; break;
      case '\"' : c = gmic_dquote; break;
      case ' ' : c = ' '; break;
      default : *(ptrd++) = '\\';
      }
      *(ptrd++) = c;
    } else if (is_dquoted) { // If non-escaped character inside string
      if (c==1) { while (c && c!=' ') c = *(++ptrs); if (!c) break; } // Discard debug info inside string
      else switch (c) {
        case '\"' : is_dquoted = false; break;
        case '$' :
          if (ptrs[1]=='?') { *(ptrd++) = '$'; is_subst = true; }
          else *(ptrd++) = gmic_dollar;
          break;
        case '{' : *(ptrd++) = gmic_lbrace; break;
        case '}' : *(ptrd++) = gmic_rbrace; break;
        case ',' : *(ptrd++) = gmic_comma; break;
        default : *(ptrd++) = c;
        }
    } else { // Non-escaped character outside string
      if (c=='\"') is_dquoted = true;
      else if (is_blank(c)) {
        *ptrd = 0;
        if (is_subst) *(++ptrd) = 1; // Item has to be substituted
        CImg<char>(item.data(),(unsigned int)(ptrd - item.data() + 1)).move_to(items);
        ptrd = item.data();
        while (is_blank(*++ptrs)) {} --ptrs; // Remove trailing spaces to next item
        is_subst = false;
      } else {
        if (c=='$' || c=='{' || c=='}' ||
            (c=='.' &&
             (ptrs==ptrs0 || is_blank(ptrs[-1]) || ptrs[-1]==',') &&
             (!ptrs[1] || is_blank(ptrs[1]) || ptrs[1]==',' ||
              (ptrs[1]=='x' && ptrs[2]>='0' && ptrs[2]<='9') ||
              (ptrs[1]=='.' &&
               (!ptrs[2] || is_blank(ptrs[2]) || ptrs[2]==',' ||
                (ptrs[2]=='x' && ptrs[3]>='0' && ptrs[3]<='9') ||
                (ptrs[2]=='.' &&
                 (!ptrs[3] || is_blank(ptrs[3]) || ptrs[3]==',' ||
                  (ptrs[3]=='x' && ptrs[4]>='0' && ptrs[4]<='9') ||
                  (ptrs[3]=='.'))))))))
          is_subst = true;
        *(ptrd++) = c;
      }
    }
  }
  if (is_dquoted) {
    CImg<char> str = CImg<char>::string(commands_line); // Discard debug info inside string
    bool _is_debug_info = false;
    const char *ptrs = ptrd = str;
    do {
      c = *(ptrs++);
      if (!c) break;
      if (c!=1) *(ptrd++) = c;
      else {
        if (!_is_debug_info) is_debug_info|=(_is_debug_info=get_debug_info(ptrs,debug_line,debug_filename));
        do { c = *(ptrs++); } while (c && c!=' ');
        if (c) ++ptrs;
      }
    } while (c);
    *ptrd = 0;
    error(true,"Invalid command line: Double quotes are not closed, in expression '%s'.",
          str.data());
  }
  if (ptrd!=item.data() && !is_blank(c)) {
    *ptrd = 0;
    if (is_subst) *(++ptrd) = 1;  // Item has to be substituted
    CImg<char>(item.data(),(unsigned int)(ptrd - item.data() + 1)).move_to(items);
  }
  if (is_debug) {
    debug("Decompose command line into %u items: ",items.size());
    cimglist_for(items,l) {
      if (items(l,0)==1) {
        if (items(l,1)) debug("  item[%u] = (debug info 0x%s)",l,items[l].data() + 1);
        else debug("  item[%u] = (undefined debug info)",l);
      } else debug("  item[%u] = '%s'",l,items[l].data());
    }
  }
  return items;
}

// Print log message.
//-------------------
gmic& gmic::print(const char *format, ...) {
  if (verbosity<1 && !is_debug) return *this;
  va_list ap;
  va_start(ap,format);
  CImg<char> message(65536);
  message[message.width() - 2] = 0;
  cimg_vsnprintf(message,message.width(),format,ap);
  strreplace_fw(message);
  if (message[message.width() - 2]) cimg::strellipsize(message,message.width() - 2);
  va_end(ap);

  // Display message.
  cimg::mutex(29);
  unsigned int &nb_carriages = cimg::output()==stdout?nb_carriages_stdout:nb_carriages_default;
  const bool is_cr = *message=='\r';
  if (is_cr) std::fputc('\r',cimg::output());
  else for (unsigned int i = 0; i<nb_carriages; ++i) std::fputc('\n',cimg::output());
  nb_carriages = 1;
  std::fprintf(cimg::output(),
               "[gmic]%s %s",
               callstack2string().data(),message.data() + (is_cr?1:0));
  std::fflush(cimg::output());
  cimg::mutex(29,0);
  return *this;
}

// Print error message, and quit interpreter.
//-------------------------------------------
gmic& gmic::error(const bool output_header, const char *const format, ...) {
  va_list ap;
  va_start(ap,format);
  CImg<char> message(1024);
  message[message.width() - 2] = 0;
  cimg_vsnprintf(message,message.width(),format,ap);
  strreplace_fw(message);
  if (message[message.width() - 2]) cimg::strellipsize(message,message.width() - 2);
  va_end(ap);

  // Display message.
  const bool is_cr = *message=='\r';
  const CImg<char> s_callstack = callstack2string();
  if (verbosity>=1 || is_debug) {
    cimg::mutex(29);
    if (is_cr) std::fputc('\r',cimg::output());
    else for (unsigned int i = 0; i<nb_carriages_default; ++i) std::fputc('\n',cimg::output());
    nb_carriages_default = 1;
    if (output_header) {
      if (is_debug_info && debug_filename<commands_files.size() && debug_line!=~0U)
        std::fprintf(cimg::output(),"[gmic]%s %s%s*** Error (file '%s', %sline #%u) *** %s%s",
                     s_callstack.data(),cimg::t_red,cimg::t_bold,
                     commands_files[debug_filename].data(),
                     is_debug_info?"":"call from ",debug_line,message.data() + (is_cr?1:0),
                     cimg::t_normal);
      else
        std::fprintf(cimg::output(),"[gmic]%s %s%s*** Error *** %s%s",
                     s_callstack.data(),cimg::t_red,cimg::t_bold,
                     message.data() + (is_cr?1:0),cimg::t_normal);
    } else
      std::fprintf(cimg::output(),"[gmic]%s %s%s%s%s",
                   s_callstack.data(),cimg::t_red,cimg::t_bold,
                   message.data() + (is_cr?1:0),cimg::t_normal);
    std::fflush(cimg::output());
    cimg::mutex(29,0);
  }

  // Store detailed error message for interpreter.
  CImg<char> full_message(512 + message.width());
  if (debug_filename<commands_files.size() && debug_line!=~0U)
    cimg_snprintf(full_message,full_message.width(),
                  "*** Error in %s (file '%s', %sline #%u) *** %s",
                  s_callstack.data(),
                  commands_files[debug_filename].data(),
                  is_debug_info?"":"call from ",debug_line,message.data() + (is_cr?1:0));
  else cimg_snprintf(full_message,full_message.width(),
                     "*** Error in %s *** %s",
                     s_callstack.data(),message.data() + (is_cr?1:0));
  CImg<char>::string(full_message).move_to(status);
  message.assign();
  is_running = false;
  throw gmic_exception(0,status);
}

// Print debug message.
//---------------------
gmic& gmic::debug(const char *format, ...) {
  if (!is_debug) return *this;
  va_list ap;
  va_start(ap,format);
  CImg<char> message(1024);
  message[message.width() - 2] = 0;
  cimg_vsnprintf(message,message.width(),format,ap);
  if (message[message.width() - 2]) cimg::strellipsize(message,message.width() - 2);
  va_end(ap);

  // Display message.
  cimg::mutex(29);
  const bool is_cr = *message=='\r';
  if (is_cr) std::fputc('\r',cimg::output());
  else for (unsigned int i = 0; i<nb_carriages_default; ++i) std::fputc('\n',cimg::output());
  nb_carriages_default = 1;

  if (is_debug_info && debug_filename<commands_files.size() && debug_line!=~0U)
    std::fprintf(cimg::output(),
                 "%s<gmic>%s#%u ",
                 cimg::t_green,callstack2string(true).data(),debug_line);
  else
    std::fprintf(cimg::output(),
                 "%s<gmic>%s ",
                 cimg::t_green,callstack2string(true).data());

  for (char *s = message.data() + (is_cr?1:0); *s; ++s) {
    char c = *s;
    if (c>=gmic_dollar && c<=gmic_dquote) switch (c) {
      case gmic_dollar : std::fprintf(cimg::output(),"\\$"); break;
      case gmic_lbrace : std::fprintf(cimg::output(),"\\{"); break;
      case gmic_rbrace : std::fprintf(cimg::output(),"\\}"); break;
      case gmic_comma : std::fprintf(cimg::output(),"\\,"); break;
      case gmic_dquote : std::fprintf(cimg::output(),"\\\""); break;
      default : std::fputc(c,cimg::output());
      } else std::fputc(c,cimg::output());
  }
  std::fprintf(cimg::output(),
               "%s",
               cimg::t_normal);
  std::fflush(cimg::output());
  cimg::mutex(29,0);
  return *this;
}

// Get variable value.
//--------------------
// May return an empty image, if requested variable is unknown.
// Returned image is a shared image, when possible.
// After return, 'varlength' contains the length of the variable content (without the ending '\0').
CImg<char> gmic::get_variable(const char *const name,
                              const unsigned int *const variables_sizes,
                              const CImgList<char> *const images_names,
                              unsigned int *const varlength) const {
  const bool
    is_global = *name=='_',
    is_thread_global = is_global && name[1]=='_';

  if (is_thread_global) cimg::mutex(30);

  // Check if variable slot exists.
  const unsigned int hash = hashcode(name,true);
  const int lmin = is_global || !variables_sizes?0:(int)variables_sizes[hash];
  CImgList<char> &vars = *variables[hash], &varnames = *variables_names[hash];
  CImg<unsigned int> &varlengths = *variables_lengths[hash];
  unsigned int ind = ~0U;
  for (int l = vars.width() - 1; l>=lmin; --l) if (!std::strcmp(varnames[l],name)) { ind = l; break; }

  // Get variable value.
  CImg<char> res;
  if (ind!=~0U) { // Regular variable name
    res.assign(vars[ind],true);
    if (varlength) *varlength = varlengths[ind];
    if (ind!=vars._width - 1) { // Modify slot position of variable to make it more accessible next time.
      unsigned int indm = (vars._width + ind)/2;
      vars[ind].swap(vars[indm]);
      varnames[ind].swap(varnames[indm]);
      cimg::swap(varlengths[ind],varlengths[indm]);
    }
  } else {
    if (images_names) { // Variable name may stand for an image index (highest index)
      const CImgList<char> &_images_names = *images_names;
      cimglist_rof(_images_names,l) if (_images_names[l] && !std::strcmp(_images_names[l],name)) { ind = l; break; }
    }
    if (ind!=~0U) {
      unsigned int tmp = std::max(1U,ind), l_tmp = 0;
      while (tmp) { ++l_tmp; tmp/=10; }
      res.assign(l_tmp + 1,1,1,1,0);
      cimg_snprintf(res,res.width(),"%u",ind);
      if (varlength) *varlength = res._width - 1;
    } else { // Variable name may stand for an environment variable
      const char *const env = std::getenv(name);
      if (env) {
        res.assign(CImg<char>::string(env,true,true),true);
        if (varlength) *varlength = res._width - 1;
      } else if (varlength) *varlength = 0;
    } // Otherwise, 'res' is empty
  }

  if (is_thread_global) cimg::mutex(30,0);
  return res;
}

// Set variable value.
//--------------------
// 'operation' can be { 0 (add new variable), '=' (replace or add), '.' (append), ',' (prepend),
//                      '+', '-', '*', '/', '%', '&', '|', '^', '<', '>' }
// If 'operation' is arithmetic (in { +,-,*,/,%,&,|,^,<,> }), 'value' must be ==0 and 'dvalue' must be defined.
// Return the new variable value.
const char *gmic::set_variable(const char *const name, const char operation,
                               const char *const value, const double dvalue,
                               const unsigned int *const variables_sizes) {
  const bool
    is_global = *name=='_',
    is_thread_global = is_global && name[1]=='_',
    is_arithmetic = operation && operation!='=' && operation!='.' && operation!=',';

  if (is_arithmetic && value)
    error(true,"Internal error: Invalid arguments to gmic::set_variable(); "
          "name='%s', operation=%d, value='%s' and dvalue='%g'.",
          name,operation,value,dvalue);

  const char *const s_operation = !is_arithmetic?0:
    operation=='+'?"+":operation=='-'?"-":operation=='*'?"*":operation=='/'?"/":operation=='%'?"%":
    operation=='&'?"&":operation=='|'?"|":operation=='^'?"^":operation=='<'?"<<":">>";
  char end;

  if (is_thread_global) cimg::mutex(30);

  // Check if variable already exists.
  const unsigned int hash = hashcode(name,true);
  const int lmin = is_global || !variables_sizes?0:(int)variables_sizes[hash];
  CImgList<char> &vars = *variables[hash], &varnames = *variables_names[hash];
  CImg<unsigned int> &varlengths = *variables_lengths[hash];
  unsigned int ind = ~0U;
  if (operation)
    for (int l = vars.width() - 1; l>=lmin; --l) if (!std::strcmp(varnames[l],name)) { ind = l; break; }

  // Create new variable slot, if needed.
  if (ind==~0U) {
    if (is_arithmetic) {
      if (is_thread_global) cimg::mutex(30,0);
      error(true,"Operator '%s=' on undefined variable '%s'.",
            s_operation,name);
    }
    ind = vars._width;
    vars.insert(1);
    CImg<char>::string(name).move_to(varnames);
    if (ind>=varlengths._width) varlengths.resize(std::max(8U,2*varlengths._width + 1),1,1,1,0);
    varlengths[ind] = 0;
  }

  // If arithmetic operation, get current variable value ('cvalue').
  double cvalue = 0;
  if (is_arithmetic) {
    if (cimg_sscanf(vars[ind],"%lf%c",&cvalue,&end)!=1) {
      if (is_thread_global) cimg::mutex(30,0);
      error(true,"Operator '%s=' on non-numerical variable '%s=%s'.",
            s_operation,name,vars[ind].data());
    }
  }

  // Store variable content.
  CImg<char> s_value;
  const unsigned int varwidth = vars[ind]._width;
  if (is_arithmetic) { // Assign with arithmetic operation
    if (varwidth<24 || varwidth>256) vars[ind].assign(24);
    cimg_snprintf(vars[ind],vars[ind].width(),"%.17g",
                  operation=='+'?cvalue + dvalue:
                  operation=='-'?cvalue - dvalue:
                  operation=='*'?cvalue*dvalue:
                  operation=='/'?cvalue/dvalue:
                  operation=='%'?cimg::mod(cvalue,dvalue):
                  operation=='&'?(double)((cimg_ulong)cvalue & (cimg_ulong)dvalue):
                  operation=='|'?(double)((cimg_ulong)cvalue | (cimg_ulong)dvalue):
                  operation=='^'?std::pow(cvalue,dvalue):
                  operation=='<'?(double)((cimg_long)cvalue << (unsigned int)dvalue):
                  (double)((cimg_long)cvalue >> (unsigned int)dvalue));
    varlengths[ind] = (unsigned int)std::strlen(vars[ind]);

  } else if ((!operation || operation=='=') && value && *value==gmic_store &&
             !std::strncmp(value + 1,"*store/",7) && value[8]) { // Assign from another image-encoded variable
    const char *const c_name = value + 8;
    const bool
      c_is_global = *c_name=='_',
      c_is_thread_global = c_is_global && c_name[1]=='_';
    if (c_is_thread_global && !is_thread_global) cimg::mutex(30);

    const unsigned int c_hash = hashcode(c_name,true);
    const int c_lmin = c_is_global || !variables_sizes?0:(int)variables_sizes[c_hash];
    CImgList<char> &c_vars = *variables[c_hash], &c_varnames = *variables_names[c_hash];
    CImg<unsigned int> &c_varlengths = *variables_lengths[c_hash];
    unsigned int c_ind = ~0U;
    for (int l = c_vars.width() - 1; l>=c_lmin; --l) if (!std::strcmp(c_varnames[l],c_name)) { c_ind = l; break; }
    if (c_ind!=~0U) {
      const unsigned int l_name = (unsigned int)std::strlen(name);
      c_vars[c_ind].get_resize(c_vars[c_ind]._width + l_name - (unsigned int)std::strlen(c_name),
                               1,1,1,0,0,1).move_to(s_value);
      cimg_snprintf(s_value,s_value._width,"%c*store/%s",gmic_store,name);
      if (c_ind!=c_vars._width - 1) { // Modify slot position of referenced image to make it more accessible next time
        unsigned int c_indm = (c_vars._width + c_ind)/2;
        c_vars[c_ind].swap(c_vars[c_indm]);
        c_varnames[c_ind].swap(c_varnames[c_indm]);
        cimg::swap(c_varlengths[c_ind],c_varlengths[c_indm]);
      }
      s_value.move_to(vars[ind]);
      varlengths[ind] = l_name + 8;
    } else {
      if (varwidth>0 && varwidth<24) *vars[ind] = 0; else vars[ind].assign(1,1,1,1,0);
      varlengths[ind] = 0;
    }

    if (c_is_thread_global && !is_thread_global) cimg::mutex(30,0);

  } else { // Append, prepend and assign
    unsigned int l_value = 0;
    if (value) { l_value = (unsigned int)std::strlen(value); s_value.assign(value,l_value + 1U,1,1,1,true); }
    else {
      s_value.assign(24);
      cimg_snprintf(s_value,s_value.width(),"%.17g",dvalue);
      l_value = (unsigned int)std::strlen(s_value);
    }

    if (operation=='.' || operation==',') { // Append and prepend
      const unsigned int varlength = varlengths[ind];
      if (!varwidth) CImg<char>(s_value._data,l_value + 1,1,1,1,true).move_to(vars[ind]);
      else if (l_value && operation=='.') { // Append
        if (varlength + l_value + 1>varwidth) { // Reallocation needed
          CImg<char> tmp(2*varwidth + l_value + 1);
          std::memcpy(tmp,vars[ind],varlength);
          tmp.move_to(vars[ind]);
        }
        std::memcpy(vars[ind]._data + varlength,s_value,l_value + 1);
      } else if (l_value && operation==',') { // Prepend
        if (varlength + l_value + 1>varwidth) { // Reallocation needed
          CImg<char> tmp(2*varwidth + l_value + 1);
          std::memcpy(tmp._data + l_value,vars[ind],varlength + 1);
          tmp.move_to(vars[ind]);
        } else std::memmove(vars[ind]._data + l_value,vars[ind]._data,varlength + 1);
        std::memcpy(vars[ind],s_value,l_value);
      }
      varlengths[ind]+=l_value;
    } else { // Assign
      if (s_value._width<=varwidth && varwidth<=8*s_value._width) // Replace in-place
        std::memcpy(vars[ind],s_value,s_value._width);
      else s_value.move_to(vars[ind]);
      varlengths[ind] = l_value;
    }
  }

  // Manage particular case of variable '_cpus': Set max number of threads for multi-threaded operators.
  if (!std::strcmp(name,"_cpus")) {
    int nb_cpus = 0;
    if (cimg_sscanf(vars[ind],"%d%c",&nb_cpus,&end)!=1 || nb_cpus<=0) {
      s_value.assign(8);
      nb_cpus = (int)cimg::nb_cpus();
      cimg_snprintf(s_value,s_value.width(),"%d",nb_cpus);
      CImg<char>::string(s_value).move_to(vars[ind]);
    }
#if cimg_use_openmp!=0
    omp_set_num_threads(nb_cpus);
#endif
  }

  // Modify slot position of modified/created variable to make it more accessible next time.
  if (ind!=vars._width - 1) {
    unsigned int indm = (vars._width + ind)/2;
    vars[ind].swap(vars[indm]);
    varnames[ind].swap(varnames[indm]);
    cimg::swap(varlengths[ind],varlengths[indm]);
  }

  if (is_thread_global) cimg::mutex(30,0);
  return vars[ind].data();
}

const char *gmic::set_variable(const char *const name, const CImg<unsigned char>& value,
                               const unsigned int *const variables_sizes) {
  if (!name || !value) return "";
  CImg<char> s_value((char*)value.data(),value.width(),value.height(),value.depth(),value.spectrum(),true);
  const bool
    is_global = *name=='_',
    is_thread_global = is_global && name[1]=='_';
  if (is_thread_global) cimg::mutex(30);
  const unsigned int hash = hashcode(name,true);
  const int lmin = is_global || !variables_sizes?0:(int)variables_sizes[hash];
  CImgList<char> &vars = *variables[hash], &varnames = *variables_names[hash];
  CImg<unsigned int> &varlengths = *variables_lengths[hash];
  unsigned int ind = ~0U;

  // Retrieve index of current definition.
  for (int l = vars.width() - 1; l>=lmin; --l) if (!std::strcmp(varnames[l],name)) { ind = l; break; }
  if (ind==~0U) { // Create new variable slot if needed
    ind = vars._width;
    vars.insert(1);
    CImg<char>::string(name).move_to(varnames);
    if (ind>=varlengths._width) varlengths.resize(std::max(8U,2*varlengths._width + 1),1,1,1,0);
    varlengths[ind] = 0;
  }
  s_value.move_to(vars[ind]); // Update variable
  varlengths[ind] = 7 + varnames[ind]._width;

  if (is_thread_global) cimg::mutex(30,0);
  return vars[ind].data();
}

// Add custom commands from a char* buffer.
//------------------------------------------
gmic& gmic::add_commands(const char *const data_commands, const char *const commands_file, const bool add_debug_info,
                         unsigned int *count_new, unsigned int *count_replaced, bool *const is_main_) {
  if (!data_commands || !*data_commands) return *this;
  cimg::mutex(23);
  CImg<char> s_body(256*1024), s_line(256*1024), s_name(257), debug_info(32);
  unsigned int line_number = 0, pos = 0;
  bool is_last_slash = false, _is_last_slash = false, is_newline = false;
  int hash = -1, l_debug_info = 0;
  char sep = 0;
  if (commands_file) {
    CImg<char>::string(commands_file).move_to(commands_files);
    CImgList<unsigned char> ltmp(commands_files.size()); // Update global variable '$_path_commands'.
    CImg<unsigned char> tmp;
    (commands_files>'x').move_to(tmp);
    tmp.resize(tmp.width() + 4,1,1,1,0,0,1);
    tmp[0] = 'G'; tmp[1] = 'M'; tmp[2] = 'Z'; tmp[3] = 0;
    tmp.unroll('y').move_to(ltmp);
    const char *const _command_files = "_path_commands";
    ltmp.get_serialize(false,(unsigned int)(9 + std::strlen(_command_files))).move_to(tmp);
    cimg_snprintf((char*)tmp.data(),tmp._height,"%c*store/%s",gmic_store,_command_files);
    set_variable(_command_files,tmp,0);
  }
  if (count_new) *count_new = 0;
  if (count_replaced) *count_replaced = 0;
  if (is_main_) *is_main_ = false;
  line_number = 1;

  for (const char *data = data_commands; *data; is_last_slash = _is_last_slash,
         line_number+=is_newline?1:0) {

    // Read new line.
    char *_line = s_line, *const line_end = s_line.end();
    while (*data!='\n' && *data && _line<line_end) *(_line++) = *(data++);
    if (_line<line_end) *_line = 0; else *(line_end - 1) = 0;
    if (*data=='\n') { is_newline = true; ++data; } else is_newline = false; // Skip next '\n'

    // Remove comments.
    _line = s_line;
    if (*_line=='#') *_line = 0; else do { // Remove comments
        if ((_line=std::strchr(_line,'#')) && is_blank(*(_line - 1))) { *--_line = 0; break; }
      } while (_line++);

    // Remove useless trailing spaces.
    char *linee = s_line.data() + std::strlen(s_line) - 1;
    while (linee>=s_line && is_blank(*linee)) --linee;
    *(linee + 1) = 0;
    char *lines = s_line; while (is_blank(*lines)) ++lines; // Remove useless leading spaces
    if (!*lines) continue; // Empty line

    // Check if last character is a '\'...
    _is_last_slash = false;
    for (_line = linee; _line>=lines && *_line=='\\'; --_line) _is_last_slash = !_is_last_slash;
    if (_is_last_slash) *(linee--) = 0; // ... and remove it if necessary
    if (!*lines) continue; // Empty line found
    *s_body = 0;
    const bool is_plus = *lines=='+';
    char
      *const nlines = lines + (is_plus?1:0),
      *const ns_name = s_name.data() + (is_plus?1:0);
    if (is_plus) *s_name = '+';
    *ns_name = sep = 0;

    if ((!is_last_slash && std::strchr(lines,':') && // Check for a command definition (or implicit '_main_')
         cimg_sscanf(nlines,"%255[a-zA-Z0-9_] %c%262143[^\n]",ns_name,&sep,s_body.data())>=2 &&
         (*nlines<'0' || *nlines>'9') && sep==':' && *s_body!='=') || ((*s_name=0), hash<0)) {
      const char *_s_body = s_body;
      if (sep==':') while (*_s_body && cimg::is_blank(*_s_body)) ++_s_body;
      CImg<char> body = CImg<char>::string(hash<0 && !*s_name?lines:_s_body);
      if (hash<0 && !*s_name) std::strcpy(s_name,"_main_");
      if (is_main_ && !std::strcmp(s_name,"_main_")) *is_main_ = true;
      hash = (int)hashcode(s_name,false);

      if (add_debug_info) { // Insert debug info code in body
        if (commands_files.width()<2)
          l_debug_info = cimg_snprintf(debug_info.data() + 1,debug_info.width() - 2,"%x",line_number);
        else
          l_debug_info = cimg_snprintf(debug_info.data() + 1,debug_info.width() - 2,"%x,%x",
                                            line_number,commands_files.width() - 1);
        if (l_debug_info>=debug_info.width() - 1) l_debug_info = debug_info.width() - 2;
        debug_info[0] = 1; debug_info[l_debug_info + 1] = ' ';
        ((CImg<char>(debug_info,l_debug_info + 2,1,1,1,true),body)>'x').move_to(body);
      }
      if (!search_sorted(s_name,commands_names[hash],commands_names[hash].size(),pos)) {
        commands_names[hash].insert(1,pos);
        commands[hash].insert(1,pos);
        commands_has_arguments[hash].insert(1,pos);
        if (count_new) ++*count_new;
      } else if (count_replaced) ++*count_replaced;
      CImg<char>::string(s_name).move_to(commands_names[hash][pos]);
      CImg<char>::vector((char)command_has_arguments(body)).
        move_to(commands_has_arguments[hash][pos]);
      body.move_to(commands[hash][pos]);

    } else { // Continuation of a previous line
      if (!is_last_slash) commands[hash][pos].back() = ' ';
      else --(commands[hash][pos]._width);
      const CImg<char> body = CImg<char>(lines,(unsigned int)(linee - lines + 2));
      commands_has_arguments[hash](pos,0) |= (char)command_has_arguments(body);
      if (add_debug_info && !is_last_slash) { // Insert code with debug info
        if (commands_files.width()<2)
          l_debug_info = cimg_snprintf(debug_info.data() + 1,debug_info.width() - 2,"%x",line_number);
        else
          l_debug_info = cimg_snprintf(debug_info.data() + 1,debug_info.width() - 2,"%x,%x",
                                       line_number,commands_files.width() - 1);
        if (l_debug_info>=debug_info.width() - 1) l_debug_info = debug_info.width() - 2;
        debug_info[0] = 1; debug_info[l_debug_info + 1] = ' ';
        ((commands[hash][pos],CImg<char>(debug_info,l_debug_info + 2,1,1,1,true),body)>'x').
          move_to(commands[hash][pos]);
      } else commands[hash][pos].append(body,'x'); // Insert code without debug info
    }
  }

  if (is_debug) {
    CImg<unsigned int> hdist(gmic_comslots);
    cimg_forX(hdist,i) hdist[i] = commands[i].size();
    const CImg<double> st = hdist.get_stats();
    cimg_snprintf(s_body,s_body.width(),
                  "Distribution of command hashes: [ %s ], min = %u, max = %u, mean = %g, std = %g.",
                  hdist.value_string().data(),(unsigned int)st[0],(unsigned int)st[1],st[2],
                  std::sqrt(st[3]));
    cimg::strellipsize(s_body,512,false);
    debug("%s",s_body.data());
  }
  cimg::mutex(23,0);
  return *this;
}

// Add commands from a file.
//---------------------------
gmic& gmic::add_commands(std::FILE *const file, const char *const commands_file, const bool add_debug_info,
                         unsigned int *count_new, unsigned int *count_replaced, bool *const is_main_) {
  if (!file) return *this;

  // Try reading it first as a .cimg file.
  try {
    CImg<char> buffer;
    buffer.load_cimg(file).unroll('x');
    buffer.resize(buffer.width() + 1,1,1,1,0);
    add_commands(buffer.data(),commands_file,add_debug_info,count_new,count_replaced,is_main_);
  } catch (...) { // If failed, read as a text file
    std::rewind(file);
    std::fseek(file,0,SEEK_END);
    const cimg_long siz = std::ftell(file);
    std::rewind(file);
    if (siz>0) {
      CImg<char> buffer((unsigned int)siz + 1);
      if (std::fread(buffer.data(),sizeof(char),siz,file)) {
        buffer[siz] = 0;
        add_commands(buffer.data(),commands_file,add_debug_info,count_new,count_replaced,is_main_);
      }
    }
  }
  return *this;
}

// Return subset indices from a selection string, as a 1-column vector.
//---------------------------------------------------------------------
CImg<unsigned int> gmic::selection2cimg(const char *const string, const unsigned int index_end,
                                        const CImgList<char>& names,
                                        const char *const command, const bool is_selection) {
#define _gmic_add_interval(i0,i1,step) \
  if (nb_intervals>=res._height) res.resize(3,std::max(8U,2*res._height),1,1,0); \
  off_intervals = 3*nb_intervals++; \
  res[off_intervals++] = (unsigned int)i0; \
  res[off_intervals++] = (unsigned int)i1; \
  res[off_intervals] = (unsigned int)step; \
  if ((unsigned int)i0<uindm) uindm = (unsigned int)i0; \
  if ((unsigned int)i1>uindM) uindM = (unsigned int)i1

#define _gmic_percent(ind) ind==0?0:ind==100?index_end - 1.:ind==50?(double)(index_end/2):(index_end - 1.)*ind/100.

  // Detect most common cases.
  CImg<unsigned int> res;
  if (string && !*string) return CImg<unsigned int>(); // Empty selection
  if (!string || (*string=='^' && !string[1])) { // Whole selection
    res.assign(1,index_end); cimg_forY(res,y) res[y] = (unsigned int)y; return res;
  } else if (*string>='0' && *string<='9' && !string[1]) { // Single positive digit
    const unsigned int ind = *string - '0';
    if (ind<index_end) return CImg<unsigned int>::vector(ind);
  } else if (*string=='-' && string[1]=='2' && string[2]==',' &&
             string[3]=='-' && string[4]=='1' && !string[5]) { // [-2,-1]
    if (index_end>=2) return CImg<unsigned int>::vector(index_end - 2,index_end - 1);
  } else if (*string=='-' && string[1]>='0' && string[2]<='9' && !string[2]) { // Single negative digit
    const unsigned int ind = index_end - string[1] + '0';
    if (ind<index_end) return CImg<unsigned int>::vector(ind);
  } else if (*string=='^') {
    if (string[1]=='-' && string[2]=='1' && !string[3]) { // ^-1
      res.assign(1,index_end - 1); cimg_forY(res,y) res[y] = (unsigned int)y; return res;
    }
    if (string[1]=='0' && !string[2]) { // ^0
      res.assign(1,index_end - 1); cimg_forY(res,y) res[y] = (unsigned int)y + 1; return res;
    }
  }

  // Parse list of intervals.
  const char
    *const stype = is_selection?"selection":"subset",
    *const ctypel = is_selection?"[":"",
    *const ctyper = is_selection?"]":"";

  const char *p = string;
  bool is_inverse = false;
  if (*p=='^') { ++p; is_inverse = true; }
  const char *const p0 = p;
  CImg<char> name;
  unsigned int nb_intervals = 0, off_intervals, uindm = ~0U, uindM = 0;
  do {
    double ind0, ind1;
    int read, istep = 1, _iind0, _iind1, iind0 = -1, iind1 = -1;

    if (p!=p0 && *p==',') ++p;
    if (cimg_sscanf(p,"%lf%n",&ind0,&read)==1) {
      p+=read;
      if (*p=='%') { ++p; ind0 = _gmic_percent(ind0); iind0 = (int)cimg::round(ind0); }
      else { _iind0 = (int)cimg::round(ind0); iind0 = _iind0<0?_iind0 + (int)index_end:_iind0; }
      if (iind0<0 || iind0>=(int)index_end) {
        if (!index_end) error(true,"Command '%s': Invalid %s '%s%s%s' (no item available).",
                              command,stype,ctypel,string,ctyper);
        error(true,"Command '%s': Invalid %s '%s%s%s' (contains index %d, not in range -%u...%u).",
              command,stype,ctypel,string,ctyper,iind0,index_end,index_end - 1);
      }
      iind1 = iind0;
      if (*p=='-') { // Sub-expression 'ind0-ind1'
        if (cimg_sscanf(++p,"%lf%n",&ind1,&read)==1) {
          p+=read;
          if (*p=='%') { ++p; ind1 = _gmic_percent(ind1); iind1 = (int)cimg::round(ind1); }
          else { _iind1 = (int)cimg::round(ind1); iind1 = _iind1<0?_iind1 + (int)index_end:_iind1; }
          if (iind1<0 || iind1>=(int)index_end)
            error(true,"Command '%s': Invalid %s '%s%s%s' (contains index %d, not in range -%u...%u).",
                  command,stype,ctypel,string,ctyper,iind1,index_end,index_end - 1);

          if (*p==':') { // Sub-expression 'ind0-ind1:step'
            if (cimg_sscanf(++p,"%d%n",&istep,&read)==1) {
              p+=read;
              if (istep<1) error(true,"Command '%s': Invalid %s '%s%s%s' (invalid step %d).",
                                 command,stype,ctypel,string,ctyper,istep);
            } else error(true,"Command '%s': Invalid %s '%s%s%s' (syntax error after colon ':').",
                         command,stype,ctypel,string,ctyper);
          }
        }
      }
      if (iind0>iind1) cimg::swap(iind0,iind1);
      _gmic_add_interval(iind0,iind1,istep);

    } else if (cimg_sscanf(p,"%255[a-zA-Z0-9_]%n",name.assign(256).data(),&read)==1 && (*name<'0' || *name>'9')) {
      p+=read;
      bool is_label_found = false;
      cimglist_for(names,l)
        if (names[l] && !std::strcmp(names[l],name)) { _gmic_add_interval(l,l,1); is_label_found = true; }
      if (!is_label_found)
        error(true,"Command '%s': Invalid %s '%s%s%s' (undefined label '%s').",
              command,stype,ctypel,string,ctyper,name.data());
    }
    if (*p && *p!=',')
      error(true,"Command '%s': Invalid %s '%s%s%s' (cannot parse '%s').",
            command,stype,ctypel,string,ctyper,p);
  } while (*p);

  // Convert list of intervals to list of indices.
  CImg<char> is_selected;
  unsigned int uind0, uind1, ustep, siz, off;

  if (res) {
    if (!is_inverse && res.height()==1) { // Single interval: optimized code
      uind0 = res[0]; uind1 = res[1]; ustep = res[2]; siz = (uind1 - uind0 + 1)/ustep;
      res.assign(1,siz); off = 0;
      for (unsigned int uind = uind0; uind<=uind1; uind+=ustep) res[off++] = uind;
    } else if (is_inverse) {
      is_selected.assign(1,index_end,1,1,1);
      for (unsigned int y = 0; y<nb_intervals; ++y) {
        uind0 = res(0,y); uind1 = res(1,y); ustep = res(2,y); siz = (uind1 - uind0 + 1)/ustep;
        for (unsigned int uind = uind0; uind<=uind1; uind+=ustep) is_selected[uind] = 0;
      }
      res.assign(1,(unsigned int)is_selected.sum()); off = 0;
      cimg_forY(is_selected,y) if (is_selected[y]) res[off++] = y;
    } else {
      is_selected.assign(1,uindM - uindm + 1,1,1,0);
      for (unsigned int y = 0; y<nb_intervals; ++y) {
        uind0 = res(0,y); uind1 = res(1,y); ustep = res(2,y); siz = (uind1 - uind0 + 1)/ustep;
        for (unsigned int uind = uind0; uind<=uind1; uind+=ustep) is_selected[uind - uindm] = 1;
      }
      res.assign(1,(unsigned int)is_selected.sum()); off = 0;
      cimg_forY(is_selected,y) if (is_selected[y]) res[off++] = y + uindm;
    }
  }
  return res;
}

// Return selection or filename strings from a set of indices.
//------------------------------------------------------------
// output_type can be { 0=display indices without brackets | 1=display indices with brackets | 2=display image names }
CImg<char>& gmic::selection2string(const CImg<unsigned int>& selection,
                                   const CImgList<char>& images_names,
                                   const unsigned int output_type,
                                   CImg<char>& res) const {
  res.assign(256);
  if (output_type<2) {
    const char *const bl = output_type?"[":"", *const br = output_type?"]":"";
    switch (selection.height()) {
    case 0:
      cimg_snprintf(res.data(),res.width()," %s%s",bl,br);
      break;
    case 1:
      cimg_snprintf(res.data(),res.width()," %s%u%s",
                    bl,selection[0],br);
      break;
    case 2:
      cimg_snprintf(res.data(),res.width(),"s %s%u,%u%s",
                    bl,selection[0],selection[1],br);
      break;
    case 3:
      cimg_snprintf(res.data(),res.width(),"s %s%u,%u,%u%s",
                    bl,selection[0],selection[1],selection[2],br);
      break;
    case 4:
      cimg_snprintf(res.data(),res.width(),"s %s%u,%u,%u,%u%s",
                    bl,selection[0],selection[1],selection[2],selection[3],br);
      break;
    case 5:
      cimg_snprintf(res.data(),res.width(),"s %s%u,%u,%u,%u,%u%s",
                    bl,selection[0],selection[1],selection[2],selection[3],selection[4],br);
      break;
    case 6:
      cimg_snprintf(res.data(),res.width(),"s %s%u,%u,%u,%u,%u,%u%s",
                    bl,selection[0],selection[1],selection[2],
                    selection[3],selection[4],selection[5],br);
      break;
    case 7:
      cimg_snprintf(res.data(),res.width(),"s %s%u,%u,%u,%u,%u,%u,%u%s",
                    bl,selection[0],selection[1],selection[2],selection[3],
                    selection[4],selection[5],selection[6],br);
      break;
    default:
      cimg_snprintf(res.data(),res.width(),"s %s%u,%u,%u,(...),%u,%u,%u%s",
                    bl,selection[0],selection[1],selection[2],
                    selection[selection.height() - 3],
                    selection[selection.height() - 2],
                    selection[selection.height() - 1],br);
    }
    return res;
  }

  switch (selection.height()) {
  case 0:
    *res = 0;
    break;
  case 1:
    cimg_snprintf(res.data(),res.width(),"%s",
                  basename(images_names[selection[0]]));
    break;
  case 2:
    cimg_snprintf(res.data(),res.width(),"%s, %s",
                  basename(images_names[selection[0]]),
                  basename(images_names[selection[1]]));
    break;
  case 3:
    cimg_snprintf(res.data(),res.width(),"%s, %s, %s",
                  basename(images_names[selection[0]]),
                  basename(images_names[selection[1]]),
                  basename(images_names[selection[2]]));
    break;
  case 4:
    cimg_snprintf(res.data(),res.width(),"%s, %s, %s, %s",
                  basename(images_names[selection[0]]),
                  basename(images_names[selection[1]]),
                  basename(images_names[selection[2]]),
                  basename(images_names[selection[3]]));
    break;
  default:
    cimg_snprintf(res.data(),res.width(),"%s, (...), %s",
                  basename(images_names[selection[0]]),
                  basename(images_names[selection.back()]));
  }
  return res;
}

// Print log message.
//-------------------
template<typename T>
gmic& gmic::print(const CImgList<T>& list, const CImg<unsigned int> *const callstack_selection,
                  const char *format, ...) {
  if (verbosity<1 && !is_debug) return *this;
  va_list ap;
  va_start(ap,format);
  CImg<char> message(65536);
  message[message.width() - 2] = 0;
  cimg_vsnprintf(message,message.width(),format,ap);
  strreplace_fw(message);
  if (message[message.width() - 2]) cimg::strellipsize(message,message.width() - 2);
  va_end(ap);

  // Display message.
  cimg::mutex(29);
  unsigned int &nb_carriages = cimg::output()==stdout?nb_carriages_stdout:nb_carriages_default;
  const bool is_cr = *message=='\r';
  if (is_cr) std::fputc('\r',cimg::output());
  else for (unsigned int i = 0; i<nb_carriages; ++i) std::fputc('\n',cimg::output());
  nb_carriages = 1;
  if (!callstack_selection || *callstack_selection)
    std::fprintf(cimg::output(),
                 "[gmic]-%u%s %s",
                 list.size(),callstack2string(callstack_selection).data(),message.data() + (is_cr?1:0));
  else std::fprintf(cimg::output(),"%s",message.data() + (is_cr?1:0));
  std::fflush(cimg::output());
  cimg::mutex(29,0);
  return *this;
}

// Print warning message.
//-----------------------
template<typename T>
gmic& gmic::warn(const CImgList<T>& list, const CImg<unsigned int> *const callstack_selection,
                 const bool force_visible, const char *const format, ...) {
  if (!force_visible && verbosity<1 && !is_debug) return *this;
  va_list ap;
  va_start(ap,format);
  CImg<char> message(1024);
  message[message.width() - 2] = 0;
  cimg_vsnprintf(message,message.width(),format,ap);
  strreplace_fw(message);
  if (message[message.width() - 2]) cimg::strellipsize(message,message.width() - 2);
  va_end(ap);

  // Display message.
  const CImg<char> s_callstack = callstack2string(callstack_selection);
  cimg::mutex(29);
  unsigned int &nb_carriages = cimg::output()==stdout?nb_carriages_stdout:nb_carriages_default;
  const bool is_cr = *message=='\r';
  if (is_cr) std::fputc('\r',cimg::output());
  else for (unsigned int i = 0; i<nb_carriages; ++i) std::fputc('\n',cimg::output());
  nb_carriages = 1;
  if (!callstack_selection || *callstack_selection) {
    if (debug_filename<commands_files.size() && debug_line!=~0U)
      std::fprintf(cimg::output(),
                   "[gmic]-%u%s %s%s*** Warning (file '%s', %sline #%u) *** %s%s",
                   list.size(),s_callstack.data(),cimg::t_magenta,cimg::t_bold,
                   commands_files[debug_filename].data(),
                   is_debug_info?"":"call from ",debug_line,message.data() + (is_cr?1:0),
                   cimg::t_normal);
    else
      std::fprintf(cimg::output(),
                   "[gmic]-%u%s %s%s*** Warning *** %s%s",
                   list.size(),s_callstack.data(),cimg::t_magenta,cimg::t_bold,
                   message.data() + (is_cr?1:0),cimg::t_normal);
  } else std::fprintf(cimg::output(),"%s%s*** Warning *** %s%s",
                      cimg::t_magenta,cimg::t_bold,message.data() + (is_cr?1:0),cimg::t_normal);
  std::fflush(cimg::output());
  cimg::mutex(29,0);
  return *this;
}

// Print error message, and quit interpreter.
//-------------------------------------------
template<typename T>
gmic& gmic::error(const bool output_header, const CImgList<T>& list,
                  const CImg<unsigned int> *const callstack_selection,
                  const char *const command, const char *const format, ...) {
  va_list ap;
  va_start(ap,format);
  CImg<char> message(1024);
  message[message.width() - 2] = 0;
  cimg_vsnprintf(message,message.width(),format,ap);
  strreplace_fw(message);
  if (message[message.width() - 2]) cimg::strellipsize(message,message.width() - 2);
  va_end(ap);

  // Display message.
  const bool is_cr = *message=='\r';
  const CImg<char> s_callstack = callstack2string(callstack_selection);
  if (verbosity>=1 || is_debug) {
    cimg::mutex(29);
    if (is_cr) std::fputc('\r',cimg::output());
    else for (unsigned int i = 0; i<nb_carriages_default; ++i) std::fputc('\n',cimg::output());
    nb_carriages_default = 1;
    if (!callstack_selection || *callstack_selection) {
      if (output_header) {
        if (debug_filename<commands_files.size() && debug_line!=~0U)
          std::fprintf(cimg::output(),
                       "[gmic]-%u%s %s%s*** Error (file '%s', %sline #%u) *** %s%s",
                       list.size(),s_callstack.data(),cimg::t_red,cimg::t_bold,
                       commands_files[debug_filename].data(),
                       is_debug_info?"":"call from ",debug_line,message.data() + (is_cr?1:0),
                       cimg::t_normal);
        else
          std::fprintf(cimg::output(),
                       "[gmic]-%u%s %s%s*** Error *** %s%s",
                       list.size(),s_callstack.data(),cimg::t_red,cimg::t_bold,
                       message.data() + (is_cr?1:0),cimg::t_normal);
      } else
        std::fprintf(cimg::output(),
                     "[gmic]-%u%s %s%s%s%s",
                     list.size(),s_callstack.data(),cimg::t_red,cimg::t_bold,
                     message.data() + (is_cr?1:0),cimg::t_normal);
    } else std::fprintf(cimg::output(),
                        "%s%s*** Error *** %s%s",
                        cimg::t_red,cimg::t_bold,
                        message.data() + (is_cr?1:0),cimg::t_normal);
    std::fflush(cimg::output());
    cimg::mutex(29,0);
  }

  // Store detailed error message for interpreter.
  CImg<char> full_message(512 + message.width());
  if (debug_filename<commands_files.size() && debug_line!=~0U)
    cimg_snprintf(full_message,full_message.width(),
                  "*** Error in %s (file '%s', %sline #%u) *** %s",
                  s_callstack.data(),
                  commands_files[debug_filename].data(),
                  is_debug_info?"":"call from ",debug_line,message.data() + (is_cr?1:0));
  else cimg_snprintf(full_message,full_message.width(),
                     "*** Error in %s *** %s",
                     s_callstack.data(),message.data() + (is_cr?1:0));
  CImg<char>::string(full_message).move_to(status);
  message.assign();
  is_running = false;
  throw gmic_exception(command,status);
}

template<typename T>
bool gmic::check_cond(const char *const expr, CImgList<T>& images, const char *const command) {
  CImg<T> &img = images.size()?images.back():CImg<T>::empty();
  bool res = false;
  float _resu = 0;
  if (!expr || !*expr) return false;
  if (img.__eval(expr,_resu)) return (bool)_resu;
  CImg<char> _expr(expr,(unsigned int)std::strlen(expr) + 1);
  strreplace_fw(_expr);
  try { if (img.eval(_expr,0,0,0,0,&images)) res = true; }
  catch (CImgException &e) {
    const char *const e_ptr = std::strstr(e.what(),": ");
    error(true,images,0,command,
          "Command '%s': Invalid argument '%s': %s",
          command,cimg::strellipsize(_expr,64,false),e_ptr?e_ptr + 2:e.what());
  }
  return res;
}

#define arg_error(command) gmic::error(true,images,0,command,"Command '%s': Invalid argument '%s'.",\
                                       command,gmic_argument_text())

// Print debug message.
//---------------------
template<typename T>
gmic& gmic::debug(const CImgList<T>& list, const char *format, ...) {
  if (!is_debug) return *this;
  va_list ap;
  va_start(ap,format);
  CImg<char> message(1024);
  message[message.width() - 2] = 0;
  cimg_vsnprintf(message,message.width(),format,ap);
  if (message[message.width() - 2]) cimg::strellipsize(message,message.width() - 2);
  va_end(ap);

  // Display message.
  cimg::mutex(29);
  const bool is_cr = *message=='\r';
  if (is_cr) std::fputc('\r',cimg::output());
  else for (unsigned int i = 0; i<nb_carriages_default; ++i) std::fputc('\n',cimg::output());
  nb_carriages_default = 1;
  if (is_debug_info && debug_filename!=~0U && debug_line!=~0U)
    std::fprintf(cimg::output(),
                 "%s<gmic>-%u%s#%u ",
                 cimg::t_green,list.size(),callstack2string(true).data(),debug_line);
  else
    std::fprintf(cimg::output(),
                 "%s<gmic>-%u%s ",
                 cimg::t_green,list.size(),callstack2string(true).data());
  for (char *s = message.data() + (is_cr?1:0); *s; ++s) {
    char c = *s;
    if (c>=gmic_dollar && c<=gmic_dquote) switch (c) {
      case gmic_dollar : std::fprintf(cimg::output(),"\\$"); break;
      case gmic_lbrace : std::fprintf(cimg::output(),"\\{"); break;
      case gmic_rbrace : std::fprintf(cimg::output(),"\\}"); break;
      case gmic_comma : std::fprintf(cimg::output(),"\\,"); break;
      case gmic_dquote : std::fprintf(cimg::output(),"\\\""); break;
      default : std::fputc(c,cimg::output());
      } else std::fputc(c,cimg::output());
  }
  std::fprintf(cimg::output(),
               "%s",
               cimg::t_normal);
  std::fflush(cimg::output());
  cimg::mutex(29,0);
  return *this;
}

// Check if a shared image of the image list is safe or not.
//----------------------------------------------------------
template<typename T>
inline bool gmic_is_valid_pointer(const T *const ptr) {
#if cimg_OS==1
  const int result = access((const char*)ptr,F_OK);
  if (result==-1 && errno==EFAULT) return false;
#elif cimg_OS==2 // #if cimg_OS==1
  return !IsBadReadPtr((void*)ptr,1);
#endif // #if cimg_OS==1
  return true;
}

template<typename T>
CImg<T>& gmic::check_image(const CImgList<T>& list, CImg<T>& img) {
  check_image(list,(const CImg<T>&)img);
  return img;
}

template<typename T>
const CImg<T>& gmic::check_image(const CImgList<T>& list, const CImg<T>& img) {
#ifdef gmic_check_image
  if (!img.is_shared() || gmic_is_valid_pointer(img.data())) return img;
  error(true,list,0,0,"Image list contains an invalid shared image (%p,%d,%d,%d,%d) "
        "(references a deallocated buffer).",
        img.data(),img.width(),img.height(),img.depth(),img.spectrum());
#else // #ifdef gmic_check_image
  cimg::unused(list);
#endif // #ifdef gmic_check_image
  return img;
}

#define gmic_check(img) check_image(images,img)

// Remove list of images in a selection.
//---------------------------------------
template<typename T>
gmic& gmic::remove_images(CImgList<T> &images, CImgList<char> &images_names,
                          const CImg<unsigned int>& selection,
                          const unsigned int start, const unsigned int end) {
  if (start==0 && end==(unsigned int)selection.height() - 1 && selection.height()==images.width()) {
    images.assign();
    images_names.assign();
  } else for (int l = (int)end; l>=(int)start; ) {
      unsigned int eind = selection[l--], ind = eind;
      while (l>=(int)start && selection[l]==ind - 1) ind = selection[l--];
      images.remove(ind,eind);
      images_names.remove(ind,eind);
    }
  return *this;
}

// This method is shared by all constructors. It initializes all the interpreter environment.
template<typename T>
gmic& gmic::_gmic(const char *const commands_line,
                  CImgList<T>& images, CImgList<char>& images_names,
                  const char *const custom_commands, const bool include_stdlib,
                  float *const p_progress, bool *const p_is_abort) {
  static bool is_first_call = true;
  cimg_exception_mode = cimg::exception_mode();
  cimg::exception_mode(0);

  // Initialize class attributes.
  cimg::mutex(22);
  if (!builtin_commands_inds) {
    builtin_commands_inds.assign(128,2,1,1,-1);
    for (unsigned int i = 0; builtin_commands_names[i]; ++i) {
      const int c = *builtin_commands_names[i];
      if (builtin_commands_inds[c]<0) builtin_commands_inds[c] = (int)i;
      builtin_commands_inds(c,1) = (int)i;
    }
  }
  if (is_first_call) {
    try { is_display_available = (bool)CImgDisplay::screen_width(); } catch (CImgDisplayException&) { }
    is_first_call = false;
  }
  cimg::mutex(22,0);

  // Initialize instance attributes.
  cimg::srand();
  setlocale(LC_NUMERIC,"C");

  delete[] commands;
  delete[] commands_names;
  delete[] commands_has_arguments;
  delete[] _variables;
  delete[] _variables_names;
  delete[] _variables_lengths;

  commands = new CImgList<char>[gmic_comslots];
  commands_names = new CImgList<char>[gmic_comslots];
  commands_has_arguments = new CImgList<char>[gmic_comslots];
  for (unsigned int l = 0; l<gmic_comslots; ++l) {
    commands_names[l].assign();
    commands[l].assign();
    commands_has_arguments[l].assign();
  }

  _variables = new CImgList<char>[gmic_varslots];
  _variables_names = new CImgList<char>[gmic_varslots];
  _variables_lengths = new CImg<unsigned int>[gmic_varslots];
  variables = new CImgList<char>*[gmic_varslots];
  variables_names = new CImgList<char>*[gmic_varslots];
  variables_lengths = new CImg<unsigned int>*[gmic_varslots];
  for (unsigned int l = 0; l<gmic_varslots; ++l) {
    _variables[l].assign();
    variables[l] = &_variables[l];
    _variables_names[l].assign();
    variables_names[l] = &_variables_names[l];
    _variables_lengths[l].assign();
    variables_lengths[l] = &_variables_lengths[l];
  }

  commands_files.assign();
  callstack.assign();
  dowhiles.assign();
  fordones.assign();
  foreachdones.assign();
  repeatdones.assign();
  light3d.assign();

  if (is_display_available) {
    display_windows.assign(gmic_winslots);
    cimg_forX(display_windows,l) display_windows[l] = new CImgDisplay;
  } else display_windows.assign();
  status.assign();

  light3d_x = light3d_y = 0;
  light3d_z = -5e8f;
  progress = p_progress?p_progress:&_progress;
  *progress = -1;

  reference_time = cimg::time();

  nb_dowhiles = nb_fordones = nb_foreachdones = nb_repeatdones = 0;
  nb_carriages_default = nb_carriages_stdout = 0;
  debug_filename = debug_line = ~0U;

  verbosity = 0;
  network_timeout = 0;

  allow_main_ = false;
  is_change = false;
  is_debug = false;
  is_running = false;
  is_start = true;
  is_return = is_quit = false;
  is_debug_info = false;
  is_abort = p_is_abort?p_is_abort:&_is_abort;
  *is_abort = false;
  is_abort_thread = false;
  is_lbrace_command = false;
  starting_commands_line = commands_line;

  // Import standard library and custom commands.
  if (include_stdlib) add_commands(gmic::decompress_stdlib().data());
  add_commands(custom_commands);

  // Set pre-defined global variables.
  CImg<char> str(16);

#if defined(_MSC_VER) || defined(WIN32) || defined(_WIN32) || defined(__WIN32__) \
 || defined(WIN64) || defined(_WIN64) || defined(__WIN64__)
  const char *s_os = "windows";
#elif defined(linux) || defined(__linux) || defined(__linux__)
  const char *s_os = "linux";
#elif defined(__OSX__) || defined(__MACOSX__) || defined(__APPLE__)
  const char *s_os = "osx";
#elif defined(BSD) || defined(__OpenBSD__) || defined(__NetBSD__) || defined(__FreeBSD__) || defined (__DragonFly__)
  const char *s_os = "bsd";
#elif defined(unix) || defined(__unix) || defined(__unix__)
  const char *s_os = "unix";
#else
  const char *s_os = "unknown";
#endif
  set_variable("_os",0,s_os);

  set_variable("_path_rc",0,gmic::path_rc());
  set_variable("_path_user",0,gmic::path_user());

  cimg_snprintf(str,str.width(),"%u",cimg::nb_cpus());
  set_variable("_cpus",0,str.data());

  set_variable("_version",0,cimg_str2(gmic_version));

#if cimg_OS==1
  cimg_snprintf(str,str.width(),"%u",(unsigned int)getpid());
#elif cimg_OS==2 // #if cimg_OS==1
  cimg_snprintf(str,str.width(),"%u",(unsigned int)_getpid());
#else // #if cimg_OS==1
  cimg_snprintf(str,str.width(),"0");
#endif // #if cimg_OS==1
  set_variable("_pid",0,str.data());

#ifdef cimg_use_vt100
  set_variable("_vt100",0,"1");
#else
  set_variable("_vt100",0,"0");
#endif // # if cimg_use_vt100

#ifdef gmic_prerelease
  set_variable("_prerelease",0,gmic_prerelease);
#else
  set_variable("_prerelease",0,"0");
#endif // #ifdef gmic_prerelease

  const char *const s_flags =
#ifdef cimg_use_board
    ",board"
#endif
#ifdef cimg_use_curl
    ",curl"
#endif
#ifdef cimg_use_fftw3
    ",fftw3"
#endif
#if cimg_display==2
    ",gdi32"
#endif
#ifdef cimg_use_jpeg
    ",jpeg"
#endif
#ifdef cimg_use_minc2
    ",minc2"
#endif
#ifdef cimg_use_magick
    ",magick"
#endif
#ifdef cimg_use_opencv
    ",opencv"
#endif
#ifdef cimg_use_openexr
    ",openexr"
#endif
#ifdef cimg_use_openmp
    ",openmp"
#endif
#ifdef gmic_is_parallel
    ",parallel"
#endif
#ifdef cimg_use_png
    ",png"
#endif
#ifdef cimg_use_tiff
    ",tiff"
#endif
#ifdef cimg_use_vt100
    ",vt100"
#endif
#if cimg_display==1
    ",x11"
#endif
#ifdef cimg_use_xrandr
    ",xrandr"
#endif
#ifdef cimg_use_xshm
    ",xshm"
#endif
#ifdef cimg_use_zlib
    ",zlib"
#endif
    "";
  set_variable("_flags",0,s_flags + 1);
  set_variable("_pixeltype",0,cimg::type<gmic_pixel_type>::string());

  // Launch G'MIC interpreter.
  const CImgList<char> items = commands_line?commands_line_to_CImgList(commands_line):CImgList<char>::empty();
  try {
    _run(items,images,images_names,true);
  } catch (gmic_exception&) {
    print(images,0,"Abort G'MIC interpreter (caught exception).\n");
    throw;
  }
  return *this;
}
bool gmic::is_display_available = false;

// Print info on selected images.
//--------------------------------
template<typename T>
gmic& gmic::print_images(const CImgList<T>& images, const CImgList<char>& images_names,
                         const CImg<unsigned int>& selection, const bool is_header) {
  if (!images || !images_names || !selection) {
    if (is_header) print(images,0,"Print image [].");
    return *this;
  }
  const bool is_verbose = verbosity>=1 || is_debug;
  CImg<char> title(256);
  if (is_header) {
    CImg<char> gmic_selection, gmic_names;
    if (is_verbose) {
      selection2string(selection,images_names,1,gmic_selection);
      selection2string(selection,images_names,2,gmic_names);
    }
    cimg::strellipsize(gmic_names,80,false);
    print(images,0,"Print image%s = '%s'.\n",
          gmic_selection.data(),gmic_names.data());
  }

  if (is_verbose) {
    cimg_forY(selection,l) {
      const unsigned int uind = selection[l];
      const CImg<T>& img = images[uind];
      const int o_verbosity = verbosity;
      const bool o_is_debug = is_debug;
      bool is_valid = true;
      verbosity = 0;
      is_debug = false;

      try { gmic_check(img); } catch (gmic_exception&) { is_valid = false; }
      is_debug = o_is_debug;
      verbosity = o_verbosity;
      cimg_snprintf(title,title.width(),"[%u] = '%s'",
                    uind,images_names[uind].data());
      cimg::strellipsize(title,80,false);
      img.gmic_print(title,is_debug,is_valid);
    }
    nb_carriages_default = 0;
  }
  return *this;
}

// Display selected images.
//-------------------------
template<typename T>
gmic& gmic::display_images(const CImgList<T>& images, const CImgList<char>& images_names,
                           const CImg<unsigned int>& selection, unsigned int *const XYZ,
                           const bool exit_on_anykey) {
  if (!images || !images_names || !selection) { print(images,0,"Display image []."); return *this; }
  const bool is_verbose = verbosity>=1 || is_debug;
  CImg<char> gmic_selection;
  if (is_verbose) selection2string(selection,images_names,1,gmic_selection);

  // Check for available display.
  if (!is_display_available) {
    cimg::unused(exit_on_anykey);
    print(images,0,"Display image%s",gmic_selection.data());
    if (is_verbose) {
      cimg::mutex(29);
      if (XYZ) std::fprintf(cimg::output(),", from point (%u,%u,%u)",XYZ[0],XYZ[1],XYZ[2]);
      std::fprintf(cimg::output()," (console output only, no display %s).\n",
                   cimg_display?"available":"support");
      std::fflush(cimg::output());
      cimg::mutex(29,0);
      print_images(images,images_names,selection,false);
    }
    return *this;
  }

  CImgList<T> visu;
  CImgList<char> t_visu;
  CImg<bool> is_valid(1,selection.height(),1,1,true);
  cimg_forY(selection,l) {
    const CImg<T>& img = images[selection[l]];
    const int o_verbosity = verbosity;
    const bool o_is_debug = is_debug;
    verbosity = 0;
    is_debug = false;
    try { gmic_check(img); } catch (gmic_exception&) { is_valid[l] = false; }
    is_debug = o_is_debug;
    verbosity = o_verbosity;
  }

  CImg<char> s_tmp;
  cimg_forY(selection,l) {
    const unsigned int uind = selection[l];
    const CImg<T>& img = images[uind];
    if (img && is_valid[l]) visu.insert(img,~0U,true);
    else visu.insert(1);
    const char *const ext = cimg::split_filename(images_names[uind]);
    const CImg<char> str = CImg<char>::string(std::strlen(ext)>6?
                                              images_names[uind].data():
                                              basename(images_names[uind]),true,true);
    s_tmp.assign(str.width() + 16);
    cimg_snprintf(s_tmp,s_tmp.width(),"[%u] %s",uind,str.data());
    s_tmp.move_to(t_visu);
  }

  CImg<char> gmic_names;
  if (visu) selection2string(selection,images_names,2,gmic_names);
  cimg::strellipsize(gmic_names,80,false);

  print(images,0,"Display image%s = '%s'",gmic_selection.data(),gmic_names.data());
  if (is_verbose) {
    cimg::mutex(29);
    if (XYZ) std::fprintf(cimg::output(),", from point (%u,%u,%u).\n",XYZ[0],XYZ[1],XYZ[2]);
    else std::fprintf(cimg::output(),".\n");
    std::fflush(cimg::output());
    nb_carriages_default = 0;
    cimg::mutex(29,0);
  }

  if (visu) {
    CImgDisplay _disp, &disp = gmic_display_window(0)?gmic_display_window(0):_disp;
    CImg<char> title(256);
    if (visu.size()==1)
      cimg_snprintf(title,title.width(),"%s (%dx%dx%dx%d)",
                    gmic_names.data(),
                    visu[0].width(),visu[0].height(),visu[0].depth(),visu[0].spectrum());
    else
      cimg_snprintf(title,title.width(),"%s (%u)",
                    gmic_names.data(),visu.size());
    cimg::strellipsize(title,80,false);
    CImg<bool> is_shared(visu.size());
    cimglist_for(visu,l) {
      is_shared[l] = visu[l].is_shared();
      visu[l]._is_shared = images[selection[l]].is_shared();
    }
    print_images(images,images_names,selection,false);
    bool is_exit = false;
    try {
      visu._gmic_display(disp,0,&t_visu,false,'x',0.5f,XYZ,exit_on_anykey,0,true,is_exit,
                         *this,visu,t_visu);
    } catch (CImgDisplayException&) {
      try { error(true,images,0,"display","Unable to display image '%s'.",gmic_names.data()); }
      catch (gmic_exception&) {}
    }
    cimglist_for(visu,l) visu[l]._is_shared = is_shared(l);
  }
  return *this;
}

// Display plots of selected images.
//----------------------------------
template<typename T>
gmic& gmic::display_plots(const CImgList<T>& images, const CImgList<char>& images_names,
                          const CImg<unsigned int>& selection,
                          const unsigned int plot_type, const unsigned int vertex_type,
                          const double xmin, const double xmax,
                          const double ymin, const double ymax,
                          const bool exit_on_anykey) {
  if (!images || !images_names || !selection) { print(images,0,"Plot image []."); return *this; }
  const bool is_verbose = verbosity>=1 || is_debug;
  CImg<char> gmic_selection;
  if (is_verbose) selection2string(selection,images_names,1,gmic_selection);

  // Check for available display.
  if (!is_display_available) {
    cimg::unused(plot_type,vertex_type,xmin,xmax,ymin,ymax,exit_on_anykey);
    print(images,0,"Plot image%s (console output only, no display %s).\n",
          gmic_selection.data(),cimg_display?"available":"support");
    print_images(images,images_names,selection,false);
    return *this;
  }

  CImgList<unsigned int> empty_indices;
  cimg_forY(selection,l) if (!gmic_check(images[selection(l)]))
    CImg<unsigned int>::vector(selection(l)).move_to(empty_indices);
  if (empty_indices && is_verbose) {
    CImg<char> eselec;
    selection2string(empty_indices>'y',images_names,1,eselec);
    warn(images,0,false,
         "Command 'plot': Image%s %s empty.",
         eselec.data(),empty_indices.size()>1?"are":"is");
  }

  CImg<char> gmic_names;
  if (is_verbose) selection2string(selection,images_names,2,gmic_names);
  print(images,0,"Plot image%s = '%s'.",
        gmic_selection.data(),gmic_names.data());

  CImgDisplay _disp, &disp = gmic_display_window(0)?gmic_display_window(0):_disp;
  bool is_first_line = false;
  cimg_forY(selection,l) {
    const unsigned int uind = selection[l];
    const CImg<T>& img = images[uind];
    if (img) {
      if (is_verbose && !is_first_line) {
        cimg::mutex(29);
        std::fputc('\n',cimg::output());
        std::fflush(cimg::output());
        cimg::mutex(29,0);
        is_first_line = true;
      }
      img.print(images_names[uind].data());
      if (!disp) disp.assign(cimg_fitscreen(CImgDisplay::screen_width()/2,CImgDisplay::screen_height()/2,1),0,0);
      img.display_graph(disp.set_title("%s (%dx%dx%dx%d)",
                                       basename(images_names[uind]),
                                       img.width(),img.height(),img.depth(),img.spectrum()),
                        plot_type,vertex_type,0,xmin,xmax,0,ymin,ymax,exit_on_anykey);
      if (is_verbose) nb_carriages_default = 0;
    }
  }
  return *this;
}

// Substitute '{}' and '$' expressions in a string.
//--------------------------------------------------
template<typename T>
CImg<char> gmic::substitute_item(const char *const source,
                                 CImgList<T>& images, CImgList<char>& images_names,
                                 CImgList<T>& parent_images, CImgList<char>& parent_images_names,
                                 const unsigned int *const variables_sizes,
                                 const CImg<unsigned int> *const command_selection,
                                 const bool is_image_expr) {
  if (!source) return CImg<char>();
  CImg<char> substituted_items(64), inbraces, substr(40), vs;
  char *ptr_sub = substituted_items.data();
  CImg<unsigned int> _ind;
  const char dot = is_image_expr?'.':0;

  for (const char *nsource = source; *nsource; )
    if (*nsource!='{' && *nsource!='$' && *nsource!=dot) {

      // If not starting with '{', '.' or '$'.
      const char *const nsource0 = nsource;
      do { ++nsource; } while (*nsource && *nsource!='{' && *nsource!='$' && *nsource!=dot);
      CImg<char>(nsource0,(unsigned int)(nsource - nsource0),1,1,1,true).
        append_string_to(substituted_items,ptr_sub);
    } else { // '{...}', '...' or '${...}' expression found
      bool is_2dollars = false, is_braces = false, is_substituted = false;
      int ind = 0, l_inbraces = 0;
      char sep = 0;
      _ind.assign();
      *substr = 0;
      if (inbraces) *inbraces = 0; else inbraces.assign(1,1,1,1,0);

      // '.' expression ('.', '..' or '...').
      if (*nsource=='.') {
        const char *str = 0;
        unsigned int p = 0, N = 0;
        if (nsource==source || *(nsource - 1)==',') {
          if (!nsource[1] || nsource[1]==',' ||
              (nsource[1]=='x' && nsource[2]>='0' && nsource[2]<='9' &&
               cimg_sscanf(nsource + 2,"%u%c",&p,&(sep=0))==1)) { str = "[-1]"; N = 1; }
          else if (nsource[1]=='.') {
            if (!nsource[2] || nsource[2]==',' ||
                (nsource[2]=='x' && nsource[3]>='0' && nsource[3]<='9' &&
                 cimg_sscanf(nsource + 3,"%u%c",&p,&(sep=0))==1)) { str = "[-2]"; N = 2; }
            else if (nsource[2]=='.') {
              if (!nsource[3] || nsource[3]==',' ||
                  (nsource[3]=='x' && nsource[4]>='0' && nsource[4]<='9' &&
                   cimg_sscanf(nsource + 4,"%u%c",&p,&(sep=0))==1)) { str = "[-3]"; N = 3; }
            }
          }
        }
        if (N) { CImg<char>::string(str,false,true).append_string_to(substituted_items,ptr_sub); nsource+=N; }
        else CImg<char>::append_string_to(*(nsource++),substituted_items,ptr_sub);
        continue;

        // '{...}' expression.
      } else if (*nsource=='{') {
        const char *const ptr_beg = nsource + 1, *ptr_end = ptr_beg;
        char delimiter = 0;
        unsigned int p = 0;
        for (p = 1; p>0 && *ptr_end; ++ptr_end) { if (*ptr_end=='{') ++p; if (*ptr_end=='}') --p; }
        if (p) {
          CImg<char>::append_string_to(*(nsource++),substituted_items,ptr_sub);
          continue;
        }
        l_inbraces = (int)(ptr_end - ptr_beg - 1);
        if (l_inbraces>0) {
          if ((*ptr_beg!='`' || ptr_beg[1]!='`') && *ptr_beg!='/' &&
              l_inbraces>2 && ptr_beg[l_inbraces - 2]==':') {
            const char del = ptr_beg[l_inbraces - 1];
            if (del==',' || del==';' || del=='/' || del=='^' || del==' ') {
              delimiter = del;
              l_inbraces-=2;
            }
          }
          inbraces.assign(ptr_beg,l_inbraces + 1).back() = 0;
          substitute_item(inbraces,images,images_names,parent_images,parent_images_names,variables_sizes,
                          command_selection,false).move_to(inbraces);
          strreplace_fw(inbraces);
        }
        nsource+=l_inbraces + 2 + (delimiter?2:0);
        if (!delimiter) delimiter = ',';
        if (!*inbraces)
          error(true,images,0,0,
                "Item substitution '{}': Empty braces.");

        // Display window features.
        if (!is_substituted && *inbraces=='*' &&
            (!inbraces[1] ||
             (inbraces[1]>='0' && inbraces[1]<='9' && !inbraces[2]) ||
             (inbraces[1]==',' && inbraces[2]) ||
             (inbraces[1]>='0' && inbraces[1]<='9' && inbraces[2]==',' && inbraces[3]))) {

          char *feature = inbraces.data() + 1;
          unsigned int wind = 0;
          if (*feature>='0' && *feature<='9') wind = (unsigned int)(*(feature++) - '0');
          char *e_feature = 0;
          if (!*feature) {
            if (!is_display_available) { *substr = '0'; substr[1] = 0; }
            else cimg_snprintf(substr,substr.width(),"%d",
                               (int)(gmic_display_window(wind) && !gmic_display_window(wind).is_closed()));
            is_substituted = true;
          } else if (*(feature++)==',') {
            do {
              is_substituted = false;
              e_feature = std::strchr(feature,',');
              if (e_feature) *e_feature = 0;
              if (!is_display_available) { *substr = '0'; substr[1] = 0; is_substituted = true; }
              else {
                CImgDisplay &disp = gmic_display_window(wind);
                bool flush_request = false;
                if (*feature=='-' &&
                    feature[1]!='w' && feature[1]!='h' && feature[1]!='d' && feature[1]!='e' &&
                    feature[1]!='u' && feature[1]!='v' && feature[1]!='n' && feature[1]!='t') {
                  flush_request = true; ++feature;
                }
                if (!*feature) { // Empty feature
                  cimg_snprintf(substr,substr.width(),"%d",disp?(disp.is_closed()?0:1):0);
                  is_substituted = true;
                } else if (*feature && !feature[1]) switch (*feature) { // Single-char features
                  case 'u' : // Screen width
                    try {
                      cimg_snprintf(substr,substr.width(),"%d",CImgDisplay::screen_width());
                    } catch (CImgDisplayException&) {
                      *substr = '0'; substr[1] = 0;
                    }
                    is_substituted = true;
                    break;
                  case 'v' : // Screen height
                    try {
                      cimg_snprintf(substr,substr.width(),"%d",CImgDisplay::screen_height());
                    } catch (CImgDisplayException&) {
                      *substr = '0'; substr[1] = 0;
                    }
                    is_substituted = true;
                    break;
                  case 'd' : // Window width
                    cimg_snprintf(substr,substr.width(),"%d",disp.window_width());
                    is_substituted = true;
                    break;
                  case 'e' : // Window height
                    cimg_snprintf(substr,substr.width(),"%d",disp.window_height());
                    is_substituted = true;
                    break;
                  case 'w' : // Display width
                    cimg_snprintf(substr,substr.width(),"%d",disp.width());
                    is_substituted = true;
                    break;
                  case 'h' : // Display height
                    cimg_snprintf(substr,substr.width(),"%d",disp.height());
                    is_substituted = true;
                    break;
                  case 'i' : // X-coordinate of window
                    cimg_snprintf(substr,substr.width(),"%d",disp.window_x());
                    is_substituted = true;
                    break;
                  case 'j' : // Y-coordinate of window
                    cimg_snprintf(substr,substr.width(),"%d",disp.window_y());
                    is_substituted = true;
                    break;
                  case 'f' : // Is fullscreen?
                    cimg_snprintf(substr,substr.width(),"%d",disp.is_fullscreen());
                    is_substituted = true;
                    break;
                  case 'n' : // Normalization type
                    cimg_snprintf(substr,substr.width(),"%d",disp.normalization());
                    is_substituted = true;
                    break;
                  case 't' : // Window title
                    cimg_snprintf(substr,substr.width(),"%s",disp.title());
                    is_substituted = true;
                    break;
                  case 'x' : // X-coordinate of mouse pointer
                    cimg_snprintf(substr,substr.width(),"%d",disp.mouse_x());
                    is_substituted = true;
                    if (flush_request) { disp._mouse_x = -1; disp._mouse_y = -1; }
                    break;
                  case 'y' : // Y-coordinate of mouse pointer
                    cimg_snprintf(substr,substr.width(),"%d",disp.mouse_y());
                    is_substituted = true;
                    if (flush_request) { disp._mouse_x = -1; disp._mouse_y = -1; }
                    break;
                  case 'b' : // State of mouse buttons
                    cimg_snprintf(substr,substr.width(),"%d",disp.button());
                    is_substituted = true;
                    if (flush_request) disp._button = 0;
                    break;
                  case 'o' : // State of mouse wheel
                    cimg_snprintf(substr,substr.width(),"%d",disp.wheel());
                    is_substituted = true;
                    if (flush_request) disp._wheel = 0;
                    break;
                  case 'c' : // Closed state of display window
                    cimg_snprintf(substr,substr.width(),"%d",(int)disp.is_closed());
                    is_substituted = true;
                    if (flush_request) disp._is_closed = false;
                    break;
                  case 'r' : // Resize event
                    cimg_snprintf(substr,substr.width(),"%d",(int)disp.is_resized());
                    is_substituted = true;
                    if (flush_request) disp._is_resized = false;
                    break;
                  case 'm' : // Move event
                    cimg_snprintf(substr,substr.width(),"%d",(int)disp.is_moved());
                    is_substituted = true;
                    if (flush_request) disp._is_moved = false;
                    break;
                  case 'k' : // Key event
                    cimg_snprintf(substr,substr.width(),"%u",disp.key());
                    is_substituted = true;
                    if (flush_request) disp._keys[0] = 0;
                    break;
                  } else if (*feature=='w' && feature[1]=='h' && !feature[2]) { // Display width*height
                  cimg_snprintf(substr,substr.width(),"%lu",
                                (unsigned long)disp.width()*disp.height());
                  is_substituted = true;
                } else if (*feature=='d' && feature[1]=='e' && !feature[2]) { // Window width*height
                  cimg_snprintf(substr,substr.width(),"%lu",
                                (unsigned long)disp.window_width()*disp.window_height());
                  is_substituted = true;
                } else if (*feature=='u' && feature[1]=='v' && !feature[2]) { // Screen width*height
                  try {
                    cimg_snprintf(substr,substr.width(),"%lu",
                                  (unsigned long)CImgDisplay::screen_width()*CImgDisplay::screen_height());
                  } catch (CImgDisplayException&) {
                    *substr = '0'; substr[1] = 0;
                  }
                  is_substituted = true;
                }
                if (*feature && !is_substituted) { // Pressed state of specified key
                  bool &ik = disp.is_key(feature);
                  cimg_snprintf(substr,substr.width(),"%d",(int)ik);
                  is_substituted = true;
                  if (flush_request) ik = false;
                }
              }
              if (is_substituted)
                CImg<char>::string(substr,false,true).append_string_to(substituted_items,ptr_sub);
              if (e_feature) {
                *e_feature = ','; feature = e_feature + 1;
                CImg<char>::append_string_to(delimiter,substituted_items,ptr_sub);
              } else feature+=std::strlen(feature);
            } while (*feature || e_feature);
            *substr = 0; is_substituted = true;
          }
        }

        // Double-backquoted string.
        if (!is_substituted && inbraces.width()>=3 && *inbraces=='`' && inbraces[1]=='`') {
          strreplace_bw(inbraces.data() + 2);
          CImg<char>(inbraces.data() + 2,inbraces.width() - 3,1,1,1,true).
            append_string_to(substituted_items,ptr_sub);
          *substr = 0; is_substituted = true;
        }

        // Escaped string.
        if (!is_substituted && inbraces.width()>=2 && *inbraces=='/') {
          const char *s = inbraces.data() + 1;
          vs.assign(inbraces.width()*4);
          const unsigned int l = strescape(s,vs);
          CImg<char>(vs,l,1,1,1,true).
            append_string_to(substituted_items,ptr_sub);
          *substr = 0; is_substituted = true;
        }

        // Sequence of character codes.
        if (!is_substituted && inbraces.width()>=3 && *inbraces=='\'' &&
            inbraces[inbraces.width() - 2]=='\'') {
          const char *s = inbraces.data() + 1;
          if (inbraces.width()>3) {
            inbraces[inbraces.width() - 2] = 0;
            cimg::strunescape(inbraces);
            for (*substr = 0; *s; ++s) {
              cimg_snprintf(substr,substr.width(),"%d%c",(int)(unsigned char)*s,delimiter);
              CImg<char>(substr.data(),(unsigned int)std::strlen(substr),1,1,1,true).
                append_string_to(substituted_items,ptr_sub);
            }
            if (*substr) --ptr_sub;
          }
          *substr = 0; is_substituted = true;
        }

        // Image feature.
        if (!is_substituted) {
          const char *feature = inbraces;
          if (l_inbraces<=2) ind = images.width() - 1; // Single-char case
          else if (cimg_sscanf(inbraces,"%d%c",&ind,&(sep=0))==2 && sep==',') {
            if (ind<0) ind+=images.width();
            if (ind<0 || ind>=images.width()) {
              if (images.width())
                error(true,images,0,0,
                      "Item substitution '{%s}': Invalid selection [%d] (not in range -%u...%u).",
                      cimg::strellipsize(inbraces,64,false),ind,images.size(),images.size() - 1);
              else
                error(true,images,0,0,
                      "Item substitution '{%s}': Invalid selection [%d] (no image data available).",
                      cimg::strellipsize(inbraces,64,false),ind);
            }
            while (*feature!=',') ++feature;
            ++feature;
          } else if (cimg_sscanf(inbraces,"%255[a-zA-Z0-9_]%c",substr.assign(256).data(),&(sep=0))==2 && sep==',') {
            selection2cimg(substr,images.size(),images_names,"Item substitution '{name,feature}'").move_to(_ind);
            if (_ind.height()!=1)
              error(true,images,0,0,
                    "Item substitution '{%s}': Invalid selection [%s], specifies multiple images.",
                    cimg::strellipsize(inbraces,64,false),substr.data());
            ind = (int)*_ind;
            while (*feature!=',') ++feature;
            ++feature;
          } else ind = images.width() - 1;

          CImg<T> &img = ind>=0?gmic_check(images[ind]):CImg<T>::empty();
          *substr = 0;
          if (!*feature)
            error(true,images,0,0,
                  "Item substitution '{%s}': Request for empty feature.",
                  cimg::strellipsize(inbraces,64,false));

          if (!feature[1]) switch (*feature) { // Single-char feature
            case 'b' : { // Image basename
              if (ind>=0 && *images_names[ind]) {
                substr.assign(std::max(substr.width(),images_names[ind].width()));
                cimg::split_filename(images_names[ind].data(),substr);
                const char *const basename = gmic::basename(substr);
                std::memmove(substr,basename,std::strlen(basename) + 1);
                strreplace_bw(substr);
              }
              is_substituted = true;
            } break;
            case 'd' : // Image depth
              cimg_snprintf(substr,substr.width(),"%d",img.depth());
              is_substituted = true;
              break;
            case 'f' : { // Image folder name
              if (ind>=0 && *images_names[ind]) {
                substr.assign(std::max(substr.width(),images_names[ind].width()));
                std::strcpy(substr,images_names[ind]);
                const char *const basename = gmic::basename(substr);
                substr[basename - substr.data()] = 0;
                strreplace_bw(substr);
              }
              is_substituted = true;
            } break;
            case 'h' : // Image height
              cimg_snprintf(substr,substr.width(),"%d",img.height());
              is_substituted = true;
              break;
            case 'n' : // Image name
              if (ind>=0 && *images_names[ind]) {
                substr.assign(std::max(substr.width(),images_names[ind].width()));
                cimg_snprintf(substr,substr.width(),"%s",images_names[ind].data());
                strreplace_bw(substr);
              }
              is_substituted = true;
              break;
            case 's' : // Number of image channels
              cimg_snprintf(substr,substr.width(),"%d",img.spectrum());
              is_substituted = true;
              break;
            case 't' : { // Ascii string from image values
              const unsigned int siz = (unsigned int)img.size();
              if (siz) {
                unsigned int strsiz = 0;
                cimg_for(img,ptr,T) if ((unsigned char)*ptr) ++strsiz; else break;
                if (strsiz) {
                  CImg<char> text(strsiz + 1), _text = text.get_shared_points(0,strsiz - 1,0,0,0);
                  _text = CImg<T>(img.data(),strsiz,1,1,1,true);
                  text.back() = 0;
                  strreplace_bw(text);
                  _text.append_string_to(substituted_items,ptr_sub);
                }
              }
              *substr = 0; is_substituted = true;
            } break;
            case 'x' : // Image extension
              if (ind>=0 && *images_names[ind]) {
                substr.assign(std::max(substr.width(),images_names[ind].width()));
                cimg_snprintf(substr,substr.width(),"%s",
                              cimg::split_filename(images_names[ind].data()));
                strreplace_bw(substr);
              }
              is_substituted = true;
              break;
            case 'w' : // Image width
              cimg_snprintf(substr,substr.width(),"%d",img.width());
              is_substituted = true;
              break;
            case '^' : { // Sequence of all pixel values
              img.value_string(delimiter).move_to(vs);
              if (vs && *vs) { --vs._width; vs.append_string_to(substituted_items,ptr_sub); }
              *substr = 0; is_substituted = true;
            } break;
            }

          const unsigned int l_feature = (unsigned int)std::strlen(feature);
          if (!is_substituted && *feature=='@') { // Subset of values
            if (l_feature>=2) {
              if (feature[1]=='^' && !feature[2]) { // All pixel values
                img.value_string(delimiter).move_to(vs);
                if (vs && *vs) { --vs._width; vs.append_string_to(substituted_items,ptr_sub); }
                *substr = 0; is_substituted = true;
              } else {
                CImg<char> subset(feature + 1,l_feature);
                subset.back() = 0;
                CImg<T> values;
                ++feature;
                CImg<char> o_status;
                status.move_to(o_status); // Save status because 'selection2cimg' may change it
                const int o_verbosity = verbosity;
                const bool o_is_debug = is_debug;
                verbosity = 0;
                is_debug = false;
                try {
                  const CImg<unsigned int> inds = selection2cimg(subset,(unsigned int)img.size(),
                                                                 CImgList<char>::empty(),"",false);
                  values.assign(1,inds.height());
                  cimg_foroff(inds,q) values[q] = img[inds[q]];
                } catch (gmic_exception &e) {
                  const char *const e_ptr = std::strstr(e.what(),": ");
                  error(true,images,0,0,
                        "Item substitution '{%s}': %s",
                        cimg::strellipsize(inbraces,64,false),e_ptr?e_ptr + 2:e.what());
                }
                is_debug = o_is_debug;
                verbosity = o_verbosity;
                o_status.move_to(status);
                cimg_foroff(values,q) {
                  cimg_snprintf(substr,substr.width(),"%.17g",(double)values[q]);
                  CImg<char>::string(substr,true,true).
                    append_string_to(substituted_items,ptr_sub);
                  *(ptr_sub - 1) = delimiter;
                }
                if (values) --ptr_sub;
              }
            }
            *substr = 0; is_substituted = true;
          }

          if (!is_substituted && l_feature==2 && *feature=='\'' && feature[1]=='\'') { // Empty string
            *substr = 0; is_substituted = true;
          }

          if (!is_substituted) { // Other mathematical expression
            const bool is_string = l_feature>=3 && *feature=='`' && inbraces[inbraces.width() - 2]=='`';
            if (is_string) { ++feature; inbraces[inbraces.width() - 2] = 0; }
            const bool is_rounded = *feature=='_';
            if (is_rounded) ++feature;

            try {
              CImg<double> output;
              img.eval(output,feature,0,0,0,0,&images);
              if (is_string) {
                vs.assign(output.height() + 1,1,1,1).fill(output).back() = 0;
                CImg<char>::string(vs,false,true).
                  append_string_to(substituted_items,ptr_sub);
              } else {
                if (output.height()>1) { // Vector-valued result
                  output.value_string(delimiter,0,is_rounded?"%g":"%.17g").move_to(vs);
                  if (vs && *vs) { --vs._width; vs.append_string_to(substituted_items,ptr_sub); }
                } else { // Scalar result
                  if (is_rounded) cimg_snprintf(substr,substr.width(),"%g",*output);
                  else cimg_snprintf(substr,substr.width(),"%.17g",*output);
                  is_substituted = true;
                }
              }
            } catch (CImgException& e) {
              const char *const e_ptr = std::strstr(e.what(),": ");
              if (is_string) inbraces[inbraces.width() - 2] = '`';
              error(true,images,0,0,
                    "Item substitution '{%s}': %s",
                    cimg::strellipsize(inbraces,64,false),e_ptr?e_ptr + 2:e.what());
            }
          }
        }
        if (is_substituted && *substr)
          CImg<char>(substr.data(),(unsigned int)std::strlen(substr),1,1,1,true).
            append_string_to(substituted_items,ptr_sub);
        continue;

        //  '${...}' and '$${...}' expressions.
      } else if (nsource[1]=='{' || (nsource[1]=='$' && nsource[2]=='{')) {
        is_2dollars = nsource[1]=='$';
        const char *const ptr_beg = nsource + 2 + (is_2dollars?1:0), *ptr_end = ptr_beg; unsigned int p = 0;
        for (p = 1; p>0 && *ptr_end; ++ptr_end) { if (*ptr_end=='{') ++p; if (*ptr_end=='}') --p; }
        if (p) {
          CImg<char>::append_string_to(*(nsource++),substituted_items,ptr_sub);
          CImg<char>::append_string_to(*(nsource++),substituted_items,ptr_sub);
          if (is_2dollars) CImg<char>::append_string_to(*(nsource++),substituted_items,ptr_sub);
          continue;
        }
        l_inbraces = (int)(ptr_end - ptr_beg - 1);
        if (l_inbraces>0) {
          inbraces.assign(ptr_beg,l_inbraces + 1).back() = 0;
          substitute_item(inbraces,images,images_names,parent_images,parent_images_names,variables_sizes,
                          command_selection,false).move_to(inbraces);
        }
        is_braces = true;
      }

      // Substitute '$?' -> String to describes the current command selection.
      if (nsource[1]=='?') {
        if (command_selection) {
          const unsigned int substr_width = (unsigned int)substr.width();
          selection2string(*command_selection,images_names,1,substr);
          CImg<char>(substr.data(),(unsigned int)std::strlen(substr),1,1,1,true).
            append_string_to(substituted_items,ptr_sub);
          substr.assign(substr_width);
          nsource+=2;
        } else CImg<char>::append_string_to(*(nsource++),substituted_items,ptr_sub);

        // Substitute '$!' -> Number of images in the list.
      } else if (nsource[1]=='!') {
        cimg_snprintf(substr,substr.width(),"%u",images.size());
        CImg<char>(substr.data(),(unsigned int)std::strlen(substr),1,1,1,true).
          append_string_to(substituted_items,ptr_sub);
        nsource+=2;

        // Substitute '$^' -> Verbosity level.
      } else if (nsource[1]=='^') {
        cimg_snprintf(substr,substr.width(),"%d",verbosity);
        CImg<char>(substr.data(),(unsigned int)std::strlen(substr),1,1,1,true).
          append_string_to(substituted_items,ptr_sub);
        nsource+=2;

        // Substitute '$|' -> Timer value.
      } else if (nsource[1]=='|') {
        cimg_snprintf(substr,substr.width(),"%.17g",(cimg::time() - reference_time)/1000.);
        CImg<char>(substr.data(),(unsigned int)std::strlen(substr),1,1,1,true).
          append_string_to(substituted_items,ptr_sub);
        nsource+=2;

        // Substitute '$/' -> Current call stack.
      } else if (nsource[1]=='/') {
        cimglist_for(callstack,i) {
          callstack[i].append_string_to(substituted_items,ptr_sub);
          *(ptr_sub - 1) = '/';
        }
        nsource+=2;

        // Substitute '$>' and '$<' -> Forward/backward index of current loop.
      } else if (nsource[1]=='>' || nsource[1]=='<') {
        if (!nb_repeatdones && !nb_dowhiles && !nb_fordones && !nb_foreachdones)
          std::strcpy(substr,"nan");
        else {
          unsigned int loop_type = 0; // 0=repeatdones, 1=dowhiles, 2=fordones, 3=foreachdones
          for (int i = (int)callstack.size() - 1; i>=0; --i) { // Find type of latest loop
            const CImg<char>& s = callstack[i];
            if (*s=='*') {
              if (s[1]=='r') { loop_type = 0; break; }
              else if (s[1]=='d') { loop_type = 1; break; }
              else if (s[1]=='f') { loop_type = s[4]!='e'?2:3; break; }
            }
          }
          switch (loop_type) {
          case 0 : { // repeat...done
            const unsigned int *const rd = repeatdones.data(0,nb_repeatdones - 1);
            if (rd[2]==~0U && nsource[1]=='<') cimg_snprintf(substr,substr.width(),"inf");
            else cimg_snprintf(substr,substr.width(),"%u",nsource[1]=='>'?rd[1]:rd[2] - 1);
          } break;
          case 1 : { // do...while
            const unsigned int *const dw = dowhiles.data(0,nb_dowhiles - 1);
            if (nsource[1]=='>') cimg_snprintf(substr,substr.width(),"%d",dw[1]);
            else std::strcpy(substr,"nan");
          } break;
          case 2 : { // for...done
            const unsigned int *const fd = fordones.data(0,nb_fordones - 1);
            if (nsource[1]=='>') cimg_snprintf(substr,substr.width(),"%d",fd[1]);
            else std::strcpy(substr,"nan");
          } break;
          case 3 : { // foreach...done
            const unsigned int *const fed = foreachdones.data(0,nb_foreachdones - 1);
            cimg_snprintf(substr,substr.width(),"%d",nsource[1]=='>'?fed[0]:fed[1] - 1);
          } break;
          }
        }
        CImg<char>(substr.data(),(unsigned int)std::strlen(substr),1,1,1,true).
          append_string_to(substituted_items,ptr_sub);
        nsource+=2;

        // Substitute '$$command' and '$${command}' -> Source of custom command.
      } else if (nsource[1]=='$' &&
                 ((is_braces &&
                   (cimg_sscanf(inbraces,"%255[a-zA-Z0-9_]",substr.assign(257).data())==1 ||
                    (*substr='+',cimg_sscanf(inbraces,"+%255[a-zA-Z0-9_]",substr.data(1)))==1) &&
                   !inbraces[std::strlen(substr)]) ||
                  (cimg_sscanf(nsource + 2,"%255[a-zA-Z0-9_]",substr.assign(257).data())==1 ||
                   (*substr='+',cimg_sscanf(nsource + 2,"+%255[a-zA-Z0-9_]",substr.data(1)))==1)) &&
                 (*substr<'0' || *substr>'9')) {
        const CImg<char>& name = is_braces?inbraces:substr;
        const unsigned int
          hash = hashcode(name,false),
          l_name = is_braces?l_inbraces + 4:(unsigned int)std::strlen(name) + 2;
        unsigned int uind = 0;
        if (search_sorted(name,commands_names[hash],commands_names[hash].size(),uind)) {
          CImg<char> sc = CImg<char>::string(commands[hash][uind],false);
          cimg_foroff(sc,p) if (sc[p]==1 && (!p || sc[p - 1]==' ')) // Discard debug info
            while (sc[p]!=' ' && p<(cimg_ulong)sc.width()) sc[p++] = 0;
          sc.discard(CImg<char>::vector(0)).autocrop(' ').unroll('x').move_to(inbraces);
          inbraces.append_string_to(substituted_items,ptr_sub);
        }
        nsource+=l_name;

        // Substitute '$name' and '${name}' -> Variable, image index or environment variable.
      } else if ((((is_braces && cimg_sscanf(inbraces,"%255[a-zA-Z0-9_]",
                                             substr.assign(256).data())==1) &&
                   !inbraces[std::strlen(substr)]) ||
                  (cimg_sscanf(nsource + 1,"%255[a-zA-Z0-9_]",substr.assign(256).data())==1)) &&
                 (*substr<'0' || *substr>'9')) {
        const CImg<char>& name = is_braces?inbraces:substr;

        unsigned int l_value = 0;
        CImg<char> value = get_variable(name,variables_sizes,&images_names,&l_value);
        const unsigned int l_name = is_braces?l_inbraces + 3:(unsigned int)std::strlen(name) + 1;
        if (value) CImg<char>(value.data(),l_value,1,1,1,true).append_string_to(substituted_items,ptr_sub);
        nsource+=l_name;

        // Substitute '${"command"}' -> Status value after command execution.
      } else if (is_braces) {
        nsource+=l_inbraces + 3;
        if (l_inbraces>0) {
          const CImgList<char>
            ncommands_line = commands_line_to_CImgList(strreplace_fw(inbraces));
          unsigned int nposition = 0;
          CImg<char>::string("*substitute").move_to(callstack);
          CImg<unsigned int> nvariables_sizes(gmic_varslots);
          cimg_forX(nvariables_sizes,l) nvariables_sizes[l] = variables[l]->size();
          const unsigned int psize = images.size();
          _run(ncommands_line,nposition,images,images_names,parent_images,parent_images_names,
               nvariables_sizes,0,inbraces,command_selection,false);
          if (images.size()!=psize)
            error(true,images,0,0,
                  "Item substitution '${\"%s\"}': Expression incorrectly changes the number of images (from %u to %u).",
                  cimg::strellipsize(inbraces,64,false),psize,images.size());
          for (unsigned int l = 0; l<gmic_varslots/2; ++l) if (variables[l]->size()>nvariables_sizes[l]) {
              if (variables_lengths[l]->_width - nvariables_sizes[l]>variables_lengths[l]->_width/2)
                variables_lengths[l]->resize(nvariables_sizes[l],1,1,1,0);
              variables_names[l]->remove(nvariables_sizes[l],variables[l]->size() - 1);
              variables[l]->remove(nvariables_sizes[l],variables[l]->size() - 1);
            }
          callstack.remove();
          is_return = false;
        }
        if (status.width()>1)
          CImg<char>(status.data(),(unsigned int)std::strlen(status),1,1,1,true).
            append_string_to(substituted_items,ptr_sub);

        // Replace '$' by itself.
      } else CImg<char>::append_string_to(*(nsource++),substituted_items,ptr_sub);
    }
  *ptr_sub = 0;
  return CImg<char>(substituted_items.data(),(unsigned int)(ptr_sub - substituted_items.data() + 1));
}

// Main parsing procedures.
//-------------------------
template<typename T>
gmic& gmic::run(const char *const commands_line, const T& pixel_type) {
  cimg::unused(pixel_type);
  CImgList<T> images;
  CImgList<char> images_names;
  return run(commands_line,images,images_names);
}

template<typename T>
gmic& gmic::run(const char *const commands_line,
                CImgList<T> &images, CImgList<char> &images_names) {
  cimg::mutex(26);
  if (is_running)
    error(true,images,0,0,
          "An instance of G'MIC interpreter %p is already running.",
          (void*)this);
  is_running = true;
  cimg::mutex(26,0);
  starting_commands_line = commands_line;
  _run(commands_line_to_CImgList(commands_line),images,images_names,true);
  is_running = false;
  return *this;
}

template<typename T>
gmic& gmic::_run(const CImgList<char>& commands_line,
                 CImgList<T> &images, CImgList<char> &images_names,
                 const bool push_new_run) {
  CImg<unsigned int> variables_sizes(gmic_varslots,1,1,1,0);
  unsigned int position = 0;
  setlocale(LC_NUMERIC,"C");
  callstack.assign(1U);
  callstack[0].assign(2,1,1,1);
  callstack[0][0] = '.';
  callstack[0][1] = 0;
  dowhiles.assign(); nb_dowhiles = 0;
  fordones.assign(); nb_fordones = 0;
  foreachdones.assign(); nb_foreachdones = 0;
  repeatdones.assign(); nb_repeatdones = 0;
  status.assign();
  nb_carriages_default = nb_carriages_stdout = 0;
  debug_filename = ~0U;
  debug_line = ~0U;
  is_change = is_debug_info = is_debug = is_quit = is_return = false;
  is_start = true;
  is_abort_thread = false;
  *progress = -1;
  cimglist_for(commands_line,l) {
    const char *it = commands_line[l].data();
    it+=*it=='-';
    if (!std::strcmp("debug",it)) { is_debug = true; break; }
  }
  return _run(commands_line,position,images,images_names,images,images_names,variables_sizes,0,0,0,push_new_run);
}

#if defined(_MSC_VER) && !defined(_WIN64)
#pragma optimize("y", off)
#endif // #if defined(_MSC_VER) && !defined(_WIN64)

template<typename T>
gmic& gmic::_run(const CImgList<char>& commands_line, unsigned int& position,
                 CImgList<T>& images, CImgList<char>& images_names,
                 CImgList<T>& parent_images, CImgList<char>& parent_images_names,
                 const unsigned int *const variables_sizes,
                 bool *const is_noarg, const char *const parent_arguments,
                 const CImg<unsigned int> *const command_selection,
                 const bool push_new_run) {
  if (*callstack.back()!='*' && (!commands_line || position>=commands_line._width)) {
    if (is_debug) debug(images,"Return from empty command '%s/'.",
                        callstack.back().data());
    return *this;
  }

  // Add/modify current run to managed list of gmic runs.
  cimg::mutex(24);
  CImgList<void*> &grl = gmic_runs();
  unsigned int ind_run = ~0U;
  CImg<void*> gr(8);
  gr[0] = (void*)this;
  gr[1] = (void*)&images;
  gr[2] = (void*)&images_names;
  gr[3] = (void*)&parent_images;
  gr[4] = (void*)&parent_images_names;
  gr[5] = (void*)variables_sizes;
  gr[6] = (void*)command_selection;
  if (!push_new_run) // Modify data for existing run
    for (int k = grl.width() - 1; k>=0; --k) {
      CImg<void*> &_gr = grl[k];
      if (_gr && _gr[0]==this) { ind_run = k; gr[7] = _gr[7]; gr.swap(_gr); break; }
    }
  if (ind_run==~0U) { ind_run = grl._width; gr[7] = get_tid(); gr.move_to(grl); } // Insert new run
  cimg::mutex(24,0);

  typedef typename cimg::superset<T,float>::type Tfloat;
  typedef typename cimg::superset<T,cimg_long>::type Tlong;
  typedef typename cimg::last<T,cimg_long>::type longT;
  const unsigned int initial_callstack_size = callstack.size(), initial_debug_line = debug_line;

  CImgList<_gmic_parallel<T> > gmic_threads;
  CImgList<unsigned int> primitives;
  CImgList<unsigned char> g_list_uc;
  CImgList<float> g_list_f;
  CImgList<char> g_list_c;
  CImgList<T> g_list;

  CImg<unsigned int> ind, ind0, ind1;
  CImg<unsigned char> g_img_uc;
  CImg<float> vertices;
  CImg<T> g_img;

  CImg<char> name, o_status, _argument_text, _argx, _argy, _argz, _argc, _title, _indices, _message, _formula, _color,
    _command(257), _s_selection(256);
  char _c0 = 0,
    *argument_text = &_c0,
    *argx = &_c0,
    *argy = &_c0,
    *argz = &_c0,
    *argc = &_c0,
    *title = &_c0,
    *indices = &_c0,
    *message = &_c0,
    *formula = &_c0,
    *color = &_c0,
    *const command = _command.data(1),
    *s_selection = _s_selection.data();
  const char *it = 0;
  *_command = '+';

// Macros below allows to allocate memory for string variables only when necessary.
#define gmic_use_var(name,siz) (name = (name!=&_c0?name:&(*_##name.assign(siz).data() = 0)))
#define gmic_use_argument_text gmic_use_var(argument_text,81)
#define gmic_use_argx gmic_use_var(argx,256)
#define gmic_use_argy gmic_use_var(argy,256)
#define gmic_use_argz gmic_use_var(argz,256)
#define gmic_use_argc gmic_use_var(argc,256)
#define gmic_use_title gmic_use_var(title,256)
#define gmic_use_indices gmic_use_var(indices,256)
#define gmic_use_message gmic_use_var(message,1024)
#define gmic_use_formula gmic_use_var(formula,4096)
#define gmic_use_color gmic_use_var(color,4096)
#define gmic_if_flr \
  if ((!_is_get && !std::strcmp("repeat",it)) || \
      (*it=='l' && ((!std::strncmp("local",it,5) && (!it[5] || it[5]=='.' || it[5]=='[')) || \
                    (!it[1] || it[1]=='.' || it[1]=='['))) || \
      (*it=='f' && ((!_is_get && !std::strcmp("for",it)) || \
                    (!std::strncmp("foreach",it,7) && (!it[7] || it[7]=='.' || it[7]=='[')))))

#define gmic_elif_flr \
  else if (!_is_get && ((*it=='}' && !it[1]) || !std::strcmp("done",it)))

  unsigned int next_debug_line = ~0U, next_debug_filename = ~0U, is_high_connectivity, uind = 0,
    boundary = 0, pattern = 0, exit_on_anykey = 0, wind = 0, interpolation = 0, hash = 0;
  char end, sep = 0, sep0 = 0, sep1 = 0, sepx = 0, sepy = 0, sepz = 0, sepc = 0, axis = 0;
  double vmin = 0, vmax = 0, value, value0, value1, nvalue, nvalue0, nvalue1;
  bool is_cond, _is_get = false, is_end_local = false, check_elif = false, run_main_ = false;
  float opacity = 0;
  int err;

  try {

    // Init interpreter environment.
    if (images.size()<images_names.size())
      images_names.remove(images.size(),images_names.size() - 1);
    else if (images.size()>images_names.size())
      images_names.insert(images.size() - images_names.size(),CImg<char>::string("[unnamed]"));

    if (is_debug) {
      if (is_start) {
        print(images,0,"Start G'MIC interpreter (debug mode, version %s).\n",cimg_str2(gmic_version));
        debug(images,"Initial command line: '%s'.",starting_commands_line);
        commands_line_to_CImgList(starting_commands_line); // Do it twice, when debug enabled
      }
      nb_carriages_default = 2;
      debug(images,"%sEnter scope '%s/'.%s",
            cimg::t_bold,callstack.back().data(),cimg::t_normal);
      is_start = false;
    }

    // Begin command line parsing.
    const int starting_verbosity = verbosity;
    if (!commands_line && is_start) { print(images,0,"Start G'MIC interpreter."); is_start = false; }

    while (position<commands_line.size() && !is_quit && !is_return) {
      const bool
        is_first_item = !position,
        was_lbrace_command = is_lbrace_command;
      *command = *s_selection = 0;
      is_lbrace_command = false;

      // Process debug info.
      if (next_debug_line!=~0U) { debug_line = next_debug_line; next_debug_line = ~0U; }
      if (next_debug_filename!=~0U) { debug_filename = next_debug_filename; next_debug_filename = ~0U; }
      while (position<commands_line.size() && *commands_line[position]==1)
        is_debug_info|=get_debug_info(commands_line[position++].data(),debug_line,debug_filename);
      if (position>=commands_line.size()) break;
      const unsigned int position_item = position; // Save initial position

      // Check consistency of the interpreter environment.
      if (images_names.size()!=images.size())
        error(true,"List of images is in an inconsistent state (%u images for %u image names). "
              "It could be caused by conccurent threads manipulating the image list at the same time.",
              images_names.size(),images.size());
      if (!callstack)
        error(true,"G'MIC encountered a fatal error (empty call stack). "
              "Please submit a bug report, at: https://github.com/GreycLab/gmic/issues");
      if (callstack.size()>=64)
        error(true,"Call stack overflow (infinite recursion?).");

      // Substitute expressions in current item.
      const char
        *const initial_item = run_main_?"_main_":commands_line[position].data(),
        *const empty_argument = "",
        *initial_argument = empty_argument;

      unsigned int position_argument = position + 1;
      while (position_argument<commands_line.size() && *(commands_line[position_argument])==1)
        is_debug_info|=get_debug_info(commands_line[position_argument++].data(),next_debug_line,next_debug_filename);
      if (position_argument<commands_line.size()) initial_argument = commands_line[position_argument];

      CImg<char> _item, _argument;
      const bool is_subst_item = (bool)commands_line[position].back();
      if (is_subst_item)
        substitute_item(initial_item,images,images_names,parent_images,parent_images_names,
                        variables_sizes,command_selection,false).move_to(_item);
      else
        CImg<char>::string(initial_item).move_to(_item);

      char *item = _item;
      const char *argument = initial_argument;
      if ((*item==',' || (*item=='{' && was_lbrace_command)) && !item[1]) { ++position; continue; }

      // Check if current item is a known command.
      const bool
        is_hyphen = *item=='-' && item[1] && item[1]!='[' && item[1]!='.' && (item[1]!='3' || item[2]!='d'),
        is_plus = *item=='+' && item[1] && item[1]!='[' && item[1]!='.' && (item[1]!='3' || item[2]!='d');
      item+=is_hyphen || is_plus?1:0;
      bool is_get = is_plus, is_specialized_get = false;

#define _gmic_eok(i) (!item[i] || item[i]=='[' || (item[i]=='.' && (!item[i + 1] || item[i + 1]=='.')))
      unsigned int hash_custom = ~0U, ind_custom = ~0U;
      const char item0 = *item, item1 = item0?item[1]:0, item2 = item1?item[2]:0;

      // Determine if specified command is a 'built-in' command (fast check when command length is 1,2 or 3).
      bool is_builtin_command = false;
      if (!item1 && (item0=='{' || item0=='}')) // Left/right braces
        is_builtin_command = true;
      else if (!is_builtin_command && item0 && _gmic_eok(1)) { switch (item0) { // Command has length 1
        case 'a': case 'b' : case 'c' : case 'd' : case 'e' : case 'f' : case 'g' : case 'h' : case 'i' : case 'j' :
        case 'k': case 'l' : case 'm' : case 'n' : case 'o' : case 'p' : case 'q' : case 'r' : case 's' : case 't' :
        case 'u': case 'v' : case 'w' : case 'x' : case 'y' : case 'z' : case '%' : case '&' : case '*' : case '+' :
        case '-': case '/' : case '<' : case '=' : case '>' : case '^' : case '|' :
          is_builtin_command = true; break;
        }
      } else if (!is_builtin_command && item0 && item1 && _gmic_eok(2)) { switch(item0) { // Command has length 2
        case '!' : is_builtin_command = item1=='='; break; // '!='
        case '<' : is_builtin_command = item1=='<' || item1=='='; break; // '<<' and '<='
        case '=' : is_builtin_command = item1=='=' || item1=='>'; break; // '==' and '=>'
        case '>' : is_builtin_command = item1=='=' || item1=='>'; break; // '>=' and '>>'
        case 'd' : is_builtin_command = item1=='o'; break; // 'do'
        case 'e' : is_builtin_command = item1=='q'; break; // 'eq'
        case 'f' : is_builtin_command = item1=='i'; break; // 'fi'
        case 'g' : is_builtin_command = item1=='e' || item1=='t'; break; // 'ge' and 'gt'
        case 'i' : is_builtin_command = item1=='f'; break; // 'if'
        case 'l' : is_builtin_command = item1=='e' || item1=='t'; break; // 'le' and 'lt'
        case 'm' : is_builtin_command = item1=='*' || item1=='/' || item1=='v'; break; // 'm*', 'm/' and 'mv'
        case 'n' : is_builtin_command = item1=='m'; break; // 'nm'
        case 'o' : is_builtin_command = item1=='r'; break; // 'or'
        case 'r' : is_builtin_command = item1=='m' || item1=='v'; break; // 'rm' and 'rv'
        case 's' : is_builtin_command = item1=='h'; break; // 'sh'
        case 'u' : is_builtin_command = item1=='m'; break; // 'um'
        case 'w' : is_builtin_command = item1>='0' && item1<='9'; break; // 'w0'..'w9'
        }
      } else if (!is_builtin_command && item0 && item1 && item2 && _gmic_eok(3)) { // Command has length 3
        is_builtin_command|= item1=='3' && item2=='d' && // '*3d','+3d','-3d', '/3d', 'j3d', 'l3d' and 'r3d'
          (item0=='*' || item0=='+' || item0=='-' || item0=='/' || item0=='j' || item0=='l' || item0=='r');
        if (!is_builtin_command) switch (item0) {
          case 'a' : is_builtin_command = (item1=='b' && item2=='s') || // 'abs', 'add' and 'and'
                                     (item2=='d' && (item1=='d' || item1=='n')); break;
          case 'b' : is_builtin_command = item1=='s' && (item2=='l' || item2=='r'); break; // 'bsl' and 'bsr'
          case 'c' : is_builtin_command = (item1=='o' && item2=='s') || // 'cos'
              (item1=='u' && item2=='t'); break; // 'cut'
          case 'd' : is_builtin_command = item1=='i' && item2=='v'; break; // 'div'
          case 'e' : is_builtin_command = (item1=='r' && item2=='f') || // 'erf'
              (item1=='x' && item2=='p'); break; // 'exp'
          case 'f' : is_builtin_command = (item1=='f' && item2=='t') || // 'fft'
              (item1=='o' && item2=='r'); break; // 'for'
          case 'l' : is_builtin_command = item1=='o' && item2=='g'; break; // 'log'
          case 'm' : is_builtin_command = (item1=='a' && (item2=='p' || item2=='x')) || // 'map' and 'max'
              (item1=='i' && item2=='n') || // 'min'
              (item1=='o' && item2=='d') || // 'mod'
              (item1=='u' && item2=='l'); break; // 'mul'
          case 'n' : is_builtin_command = (item1=='e' && item2=='q') || // 'neq'
              (item1=='m' && item2=='d'); break; // 'nmd'
          case 'p' : is_builtin_command = item1=='o' && item2=='w'; break; // 'pow'
          case 'r' : is_builtin_command = item1=='o' && (item2=='l' || item2=='r'); break; // 'rol' and 'ror'
          case 's' : is_builtin_command = (item1=='e' && item2=='t') || // 'set'
              (item1=='i' && item2=='n') || // 'sin'
              (item1=='q' && item2=='r') || // 'sqr'
              (item1=='u' && item2=='b') || // 'sub'
              (item1=='v' && item2=='d'); break; // 'svd'
          case 't' : is_builtin_command = item1=='a' && item2=='n'; break; // 'tan'
          case 'x' : is_builtin_command = item1=='o' && item2=='r'; break; // 'xor'
          }
      }

      bool is_command = is_builtin_command;

      if (!is_builtin_command) {
        *command = sep0 = sep1 = sep = 0;

        // Extract command name.
        // (same as but faster than 'err = cimg_sscanf(item,"%255[a-zA-Z_0-9]%c%c",command,&sep0,&sep1,&sep);').
        const char *ps = item;
        char *pd = command;
        char *const pde = _command.end() - 1;
        for (err = 0; *ps && pd<pde; ++ps) {
          const char c = *ps;
          if ((c>='a' && c<='z') || (c>='A' && c<='Z') || (c>='0' && c<='9') || c=='_') *(pd++) = c;
          else break;
        }
        if (pd!=command) {
          *pd = 0;
          ++err;
          if (*ps) {
            sep0 = *(ps++);
            ++err;
            if (*ps) {
              sep1 = *(ps++);
              ++err;
              if (*ps) { sep = *(ps++); ++err; }
            }
          }
        }

        is_command =
          (err==1 || (err==2 && sep0=='.') || (err>=3 && (sep0=='[' || (sep0=='.' && sep1=='.' && sep!='=')))) &&
          (*item<'0' || *item>'9');

        if (is_command) {
          if (!is_builtin_command) { // Search for known built-in command name
            const int
              _ind0 = builtin_commands_inds[(unsigned int)*command],
              _ind1 = builtin_commands_inds((unsigned int)*command,1U);
            if (_ind0>=0)
              is_builtin_command = search_sorted(command,builtin_commands_names + _ind0,
                                                 _ind1 - _ind0 + 1U,uind);
          }
          if (!is_builtin_command) { // Search for a custom command name
            bool found_command = false;

            if (is_get) { // Search for specialization '+command' that has priority over 'command'
              hash_custom = hashcode(_command,false);
              found_command = search_sorted(_command,commands_names[hash_custom],
                                            commands_names[hash_custom].size(),ind_custom);
              if (found_command) { is_get = false; is_specialized_get = true; }
            }
            if (!found_command) { // Search for regular command
              hash_custom = hashcode(command,false);
              found_command = search_sorted(command,commands_names[hash_custom],
                                            commands_names[hash_custom].size(),ind_custom);
            }
            if (!found_command && !is_get) { // Finally, search for specialization '+command' invoked as 'command'
              hash_custom = hashcode(_command,false);
              found_command = search_sorted(_command,commands_names[hash_custom],
                                            commands_names[hash_custom].size(),ind_custom);
              if (found_command) is_specialized_get = true;
            }
            if (!found_command) hash_custom = ind_custom = ~0U;
          }
        }
        is_command|=is_builtin_command;
      }

      // Split command/selection, if necessary.
      bool is_selection = false;
      const unsigned int siz = images._width, selsiz = _s_selection._width;
      CImg<unsigned int> selection;
      if (is_command) {
        sep0 = sep1 = 0;
        strreplace_fw(item);

        // Extract selection.
        // (same as but faster than 'err = cimg_sscanf(item,"%255[^[]%c%255[a-zA-Z_0-9.eE%^,:+-]%c%c",
        //                                             command,&sep0,s_selection,&sep1,&end);
        if (selsiz<_item._width) { // Expand size for getting a possibly large selection
          _s_selection.assign(_item.width());
          s_selection = _s_selection.data();
          *s_selection = 0;
        }

        const char *ps = item;
        char *pd = command;
        char *const pde = _command.end() - 1;
        for (err = 0; *ps && *ps!='[' && pd<pde; ++ps) *(pd++) = *ps;
        if (pd!=command) {
          *pd = 0;
          ++err;
          if (*ps) {
            sep0 = *(ps++);
            ++err;
            if (*ps) {
              const char *pde2 = _s_selection.end() - 1;
              for (pd = s_selection; *ps && pd<pde2; ++ps) {
                const char c = *ps;
                if ((c>='a' && c<='z') || (c>='A' && c<='Z') || (c>='0' && c<='9') ||
                    c=='_' || c=='.' || c=='e' || c=='E' || c=='%' || c=='^' || c==','
                    || c==':' || c=='+' || c=='-') *(pd++) = c;
                else break;
              }
              if (pd!=s_selection) {
                *pd = 0;
                ++err;
                if (*ps) {
                  sep1 = *(ps++);
                  ++err;
                  if (*ps) ++err;
                }
              }
            }
          }
        }

        const unsigned int l_command = err==1?(unsigned int)std::strlen(command):0;
        if (err==1 && l_command>=2 && command[l_command - 1]=='.') { // Selection shortcut
          err = 4; sep0 = '['; sep1 = ']';
          if (command[l_command - 2]!='.') {
            *s_selection = '-'; s_selection[1] = '1'; s_selection[2] = 0; command[l_command - 1] = 0;
          }
          else if (l_command>=3 && command[l_command - 3]!='.') {
            *s_selection = '-'; s_selection[1] = '2'; s_selection[2] = 0; command[l_command - 2] = 0;
          }
          else if (l_command>=4 && command[l_command - 4]!='.') {
            *s_selection = '-'; s_selection[1] = '3'; s_selection[2] = 0; command[l_command - 3] = 0;
          }
          else { is_command = false; ind_custom = ~0U; *s_selection = 0; }
        }

        if (err==1) { // No selection -> all images
          if (!std::strcmp(command,"pass")) selection.assign(1,parent_images.size());
          else if (!((*command=='=' && command[1]=='>' && !command[2]) ||
                     (*command=='n' && command[1]=='m' && !command[2]) ||
                     !std::strcmp(command,"name"))) selection.assign(1,siz);
          cimg_forY(selection,y) selection[y] = (unsigned int)y;
        } else if (err==2 && sep0=='[' && item[std::strlen(command) + 1]==']') { // Empty selection
          is_selection = true;
        } else if (err==4 && sep1==']') { // Other selections
          is_selection = true;
          if (!is_get && (!std::strcmp("wait",command) ||
                          !std::strcmp("cursor",command)))
            selection2cimg(s_selection,10,CImgList<char>::empty(),command).move_to(selection);
          else if (*command=='i' && (!command[1] || !std::strcmp("input",command)))
            selection2cimg(s_selection,siz + 1,images_names,command,true).move_to(selection);
          else if (!is_get &&
                   ((*command=='e' && (!command[1] ||
                                       !std::strcmp("echo",command) ||
                                       !std::strcmp("error",command))) ||
                    !std::strcmp("warn",command)))
            selection2cimg(s_selection,callstack.size(),CImgList<char>::empty(),command).move_to(selection);
          else if (!std::strcmp("pass",command))
            selection2cimg(s_selection,parent_images.size(),parent_images_names,command).move_to(selection);
          else
            selection2cimg(s_selection,siz,images_names,command).move_to(selection);
        } else {
          std::strncpy(command,item,_command.width() - 2);
          is_command = false; ind_custom = ~0U;
          command[_command.width() - 2] = *s_selection = 0;
        }
      } else {
        std::strncpy(command,item,_command.width() - 2);
        command[_command.width() - 2] = *s_selection = 0;
      }
      position = position_argument;
      if (_s_selection._width!=selsiz) { // Go back to initial size for selection image.
        _s_selection.assign(selsiz);
        s_selection = _s_selection.data();
        *s_selection = 0;
      }

      const bool
        is_command_verbose = is_get?false:
          is_command && *item=='v' && (!item[1] || !std::strcmp(item,"verbose")),
        is_command_echo = is_command_verbose?false:
          is_command && *command=='e' && (!command[1] || !std::strcmp(command,"echo")),
        is_command_error = is_get || is_command_verbose || is_command_echo?false:
          is_command && *command=='e' && !std::strcmp(command,"error"),
        is_command_warn = is_command_verbose || is_command_echo || is_command_error?false:
          is_command && *command=='w' && !std::strcmp(command,"warn"),
        is_command_input = is_command_verbose || is_command_echo || is_command_error ||
          is_command_warn?false:
          is_command && *command=='i' && (!command[1] || !std::strcmp(command,"input")),
        is_command_check = is_get || is_command_verbose || is_command_echo || is_command_error ||
          is_command_warn || is_command_input?false:
          is_command && *command=='c' && !std::strcmp(item,"check"),
        is_command_skip = is_get || is_command_verbose || is_command_echo || is_command_error ||
          is_command_warn || is_command_input || is_command_check?false:
        is_command && *command=='s' && !std::strcmp(item,"skip");
      bool is_subst_arg = position_argument<commands_line.size()?(bool)commands_line[position_argument].back():false;

      // Check for verbosity command, prior to the first output of a log message.
      bool is_verbose = verbosity>=1 || is_debug, is_verbose_argument = false;
      const int prev_verbosity = verbosity;
      if (is_command_verbose) {
        // Do a first fast check.
        if (*argument=='-' && !argument[1]) { --verbosity; is_verbose_argument = true; }
        else if (*argument=='+' && !argument[1]) { ++verbosity; is_verbose_argument = true; }
        else {
          gmic_substitute_args(false);
          if (*argument=='-' && !argument[1]) { --verbosity; is_verbose_argument = true; }
          else if (*argument=='+' && !argument[1]) { ++verbosity; is_verbose_argument = true; }
          else {
            float level = 0;
            if (cimg_sscanf(argument,"%f%c",&level,&end)==1) {
              verbosity = (int)cimg::round(level);
              is_verbose_argument = true;
            }
            else arg_error("verbose");
          }
        }
      }
      is_verbose = verbosity>=1 || is_debug;
      const bool is_very_verbose = verbosity>1 || is_debug;

      // Generate strings for displaying image selections when verbosity>=1.
      CImg<char> gmic_selection;
      if (is_debug ||
          (verbosity>=1 && !is_command_check && !is_command_skip &&
           !is_command_verbose && !is_command_echo && !is_command_error && !is_command_warn))
        selection2string(selection,images_names,1,gmic_selection);

      if (is_debug) {
        if (std::strcmp(item,initial_item))
          debug(images,"Item[%u]: '%s' -> '%s', selection%s.",
                position_item,initial_item,item,gmic_selection.data());
        else
          debug(images,"Item[%u]: '%s', selection%s.",
                position_item,initial_item,gmic_selection.data());
      }

      // Display starting message.
      if (is_start) {
        print(images,0,"Start G'MIC interpreter.");
        is_start = false;
      }

      // Cancellation point.
      if (*is_abort || is_abort_thread)
        throw CImgAbortException();

      // Begin command interpretation.
      if (is_command) {

        // Convert command shortcuts to full names.
        char command0 = *command;
        const char
          command1 = command0?command[1]:0, command2 = command1?command[2]:0, command3 = command2?command[3]:0;

        static const char* onechar_shortcuts[] = {
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 0-31
          0,0,0,0,0,"mod","and",0,0,0,"mul","add",0,"sub",0,"div",0,0,0,0,0,0,0,0,0,0,0,0, // 32-59
          "lt","set","gt",0, // 60-63
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"pow",0, // 64-95
          0,"append","blur","cut","display","echo","fill",0,0,"input","image","keep", // 96-107
          "local","command","normalize","output","print","quit","resize","split","text","status", // 108-117
          "verbose","window","exec","unroll","crop",0,"or","done",0,0 // 118-127
        };

        if (!command1) { // Single-char shortcut
          const bool
            is_mquvx = command0=='m' || command0=='q' || command0=='u' || command0=='v' || command0=='x' ||
                       command0=='}',
            is_deiowx = command0=='d' || command0=='e' || command0=='i' || command0=='o' || command0=='w' ||
                        command0=='x';
          if ((unsigned int)command0<128 && onechar_shortcuts[(unsigned int)command0] &&
              !(is_mquvx && (is_get || is_selection)) && !(is_deiowx && is_get)) {
            if (is_mquvx) {
              CImg<char>::string(onechar_shortcuts[(unsigned int)command0]).move_to(_item);
              *command = 0;
            } else std::strcpy(command,onechar_shortcuts[(unsigned int)command0]);
          }

        } else if (!command2) { // Two-chars shortcuts
          if (command0=='s' && command1=='h') std::strcpy(command,"shared");
          else if (command0=='m' && command1=='v') std::strcpy(command,"move");
          else if (!is_get && ((command0=='=' && command1=='>') || (command0=='n' && command1=='m')))
            std::strcpy(command,"name");
          else if (command0=='r' && command1=='m') std::strcpy(command,"remove");
          else if (command0=='r' && command1=='v') std::strcpy(command,"reverse");
          else if (command0=='<' && command1=='<') std::strcpy(command,"bsl");
          else if (command0=='>' && command1=='>') std::strcpy(command,"bsr");
          else if (command0=='=' && command1=='=') std::strcpy(command,"eq");
          else if (command0=='>' && command1=='=') std::strcpy(command,"ge");
          else if (command0=='<' && command1=='=') std::strcpy(command,"le");
          else if (command0=='m' && command1=='/') std::strcpy(command,"mdiv");
          else if (command0=='m' && command1=='*') std::strcpy(command,"mmul");
          else if (command0=='!' && command1=='=') std::strcpy(command,"neq");
          else if (command0=='u' && command1=='m') CImg<char>::string("uncommand").move_to(_item);

        } else if (!command3 && command1=='3' && command2=='d') switch (command0) {
            // Three-chars shortcuts, ending with '3d'.
          case 'j' : std::strcpy(command,"object3d"); break;
          case '+' : std::strcpy(command,"add3d"); break;
          case '/' : std::strcpy(command,"div3d"); break;
          case 'l' : if (!is_get && !is_selection) CImg<char>::string("light3d").move_to(_item); break;
          case '*' : std::strcpy(command,"mul3d"); break;
          case 'r' : std::strcpy(command,"rotate3d"); break;
          case '-' : std::strcpy(command,"sub3d"); break;
          } else if (!is_get && !command3 && command0=='n' && command1=='m' && command2=='d') {
          std::strcpy(command,"named"); // Shortcut 'nmd' for 'named".
        }
        if (item!=_item.data() + (is_hyphen || is_plus?1:0)) item = _item;
        command0 = *command?*command:*item;

        // Dispatch to dedicated parsing code, regarding the first character of the command.
        // We rely on the compiler to optimize this using an associative array (verified with g++).
        if (!is_builtin_command) goto gmic_commands_others;
        switch (command0) {
        case 'a' : goto gmic_commands_a;
        case 'b' : goto gmic_commands_b;
        case 'c' : goto gmic_commands_c;
        case 'd' : goto gmic_commands_d;
        case 'e' : goto gmic_commands_e;
        case 'f' : goto gmic_commands_f;
        case 'g' : goto gmic_commands_g;
        case 'h' : goto gmic_commands_h;
        case 'i' : goto gmic_commands_i;
        case 'k' : goto gmic_commands_k;
        case 'l' : goto gmic_commands_l;
        case 'm' : goto gmic_commands_m;
        case 'n' : goto gmic_commands_n;
        case 'o' : goto gmic_commands_o;
        case 'p' : goto gmic_commands_p;
        case 'q' : goto gmic_commands_q;
        case 'r' : goto gmic_commands_r;
        case 's' : goto gmic_commands_s;
        case 't' : goto gmic_commands_t;
        case 'u' : goto gmic_commands_u;
        case 'v' : goto gmic_commands_v;
        case 'w' : goto gmic_commands_w;
        case 'x' : goto gmic_commands_x;
        default : goto gmic_commands_others;
        }

        //-----------------------------
        // Commands starting by 'a...'
        //-----------------------------
      gmic_commands_a :

        // Append.
        if (!std::strcmp("append",command)) {
          gmic_substitute_args(true);
          float align = 0;
          axis = sep = 0;
          if ((cimg_sscanf(argument,"%c%c",
                           &axis,&end)==1 ||
               cimg_sscanf(argument,"%c,%f%c",
                           &axis,&align,&end)==2) &&
              is_xyzc(axis)) {
            print(images,0,"Append image%s along the '%c'-axis, with alignment %g.",
                  gmic_selection.data(),
                  axis,align);
            if (selection) {
              cimg_forY(selection,l) if (gmic_check(images[selection[l]]))
                g_list.insert(gmic_check(images[selection[l]]),~0U,true);
              CImg<T> img = g_list.get_append(axis,align);
              if (is_get) {
                img.move_to(images);
                images_names[selection[0]].get_copymark().move_to(images_names);
              } else if (selection.height()>=2) {
                remove_images(images,images_names,selection,1,selection.height() - 1);
                img.move_to(images[selection[0]].assign());
              }
              g_list.assign();
            }
          } else if ((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c,%c%c",
                                  &(*gmic_use_indices=0),&sep,&axis,&end)==3 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c,%c,%f%c",
                                  indices,&sep,&axis,&(align=0),&end)==4) &&
                     is_xyzc(axis) && sep==']' &&
                     (ind=selection2cimg(indices,images.size(),images_names,"append")).height()==1) {
            print(images,0,"Append image [%u] to image%s, along the '%c'-axis, with alignment %g.",
                  *ind,gmic_selection.data(),axis,align);
            const CImg<T> img0 = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply(append(img0,axis,align),false);
          } else arg_error("append");
          is_change = true;
          ++position;
          continue;
        }

        // Autocrop.
        if (!std::strcmp("autocrop",command)) {
          gmic_substitute_args(false);
          if (*argument && cimg_sscanf(argument,"%4095[0-9.,eEinfa+-]%c",gmic_use_formula,&end)==1)
            try { CImg<T>(1).fill_from_values(argument,true).move_to(g_img); }
            catch (CImgException&) { g_img.assign(); }
          if (g_img) {
            print(images,0,"Auto-crop image%s by vector '%s'.",
                  gmic_selection.data(),
                  gmic_argument_text_printed());
            ++position;
          } else print(images,0,"Auto-crop image%s.",
                       gmic_selection.data());
          cimg_forY(selection,l) {
            if (g_img) {
              CImg<T>& img = images[selection[l]];
              g_img.assign(img.spectrum()).fill_from_values(argument,true);
              gmic_apply(gmic_autocrop(g_img),false);
            } else gmic_apply(gmic_autocrop(),false);
          }
          g_img.assign();
          is_change = true;
          continue;
        }

        // Add.
        gmic_arithmetic_command("add",
                                operator+=,
                                "Add %g%s to image%s",
                                value,ssep,gmic_selection.data(),Tfloat,
                                operator+=,
                                "Add image [%d] to image%s",
                                ind[0],gmic_selection.data(),
                                operator_pluseq,
                                "Add expression %s to image%s",
                                gmic_argument_text_printed(),gmic_selection.data(),
                                "Add image%s");

        // Add 3D objects together, or shift a 3D object.
        if (!std::strcmp("add3d",command)) {
          gmic_substitute_args(true);
          double tx = 0, ty = 0, tz = 0;
          sep = *indices = 0;
          if (cimg_sscanf(argument,"%lf%c",
                          &tx,&end)==1 ||
              cimg_sscanf(argument,"%lf,%lf%c",
                          &tx,&ty,&end)==2 ||
              cimg_sscanf(argument,"%lf,%lf,%lf%c",
                          &tx,&ty,&tz,&end)==3) {
            print(images,0,"Shift 3D object%s by displacement (%g,%g,%g).",
                  gmic_selection.data(),
                  tx,ty,tz);
            cimg_forY(selection,l) {
              uind = selection[l];
              CImg<T>& img = images[uind];
              try { gmic_apply(shift_CImg3d((float)tx,(float)ty,(float)tz),true); }
              catch (CImgException&) {
                if (!img.is_CImg3d(true,&(*gmic_use_message=0)))
                  error(true,images,0,0,
                        "Command 'add3d': Invalid 3D object [%d], in image%s (%s).",
                        uind,gmic_selection_err.data(),message);
                else throw;
              }
            }
            ++position;
          } else if (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep,&end)==2 &&
                     sep==']' &&
                     (ind=selection2cimg(indices,images.size(),images_names,"add3d")).height()==1) {
            const CImg<T> img0 = gmic_image_arg(*ind);
            print(images,0,"Merge 3D object%s with 3D object [%u].",
                  gmic_selection.data(),*ind);
            cimg_forY(selection,l) {
              const unsigned int _ind = selection[l];
              CImg<T>& img = gmic_check(images[_ind]);
              g_list.assign(2);
              g_list[0].assign(img,true);
              g_list[1].assign(img0,true);
              CImg<T> res;
              try { CImg<T>::append_CImg3d(g_list).move_to(res); }
              catch (CImgException&) {
                if (!img0.is_CImg3d(true,&(*gmic_use_message=0)))
                  error(true,images,0,0,
                        "Command 'add3d': Invalid 3D object [%u], in specified "
                        "argument '%s' (%s).",
                        *ind,gmic_argument_text(),message);
                else if (!img.is_CImg3d(true,message))
                  error(true,images,0,0,
                        "Command 'add3d': Invalid 3D object [%d], in image%s (%s).",
                        _ind,gmic_selection_err.data(),message);
                else throw;
              }
              if (is_get) {
                res.move_to(images);
                images_names[_ind].get_copymark().move_to(images_names);
              } else res.move_to(images[_ind].assign());
            }
            g_list.assign();
            ++position;
          } else {
            print(images,0,"Merge 3D object%s.",
                  gmic_selection.data());
            if (selection) {
              g_list.assign(selection.height());
              cimg_forY(selection,l) g_list[l].assign(gmic_check(images[selection[l]]),true);
              CImg<T> img;
              try { CImg<T>::append_CImg3d(g_list).move_to(img); }
              catch (CImgException&) {
                cimg_forY(selection,l) {
                  uind = selection[l];
                  if (!images[uind].is_CImg3d(true,&(*gmic_use_message=0)))
                    error(true,images,0,0,
                          "Command 'add3d': Invalid 3D object [%d], in image%s (%s).",
                          uind,gmic_selection_err.data(),message);
                }
                throw;
              }
              if (is_get) {
                img.move_to(images);
                images_names[selection[0]].get_copymark().move_to(images_names);
              } else if (selection.height()>=2) {
                remove_images(images,images_names,selection,1,selection.height() - 1);
                img.move_to(images[selection[0]].assign());
              }
              g_list.assign();
            }
          }
          is_change = true;
          continue;
        }

        // Absolute value.
        gmic_simple_command("abs",abs,"Compute pointwise absolute value of image%s.");

        // Bitwise and.
        gmic_arithmetic_command("and",
                                operator&=,
                                "Compute bitwise AND of image%s by %g%s",
                                gmic_selection.data(),value,ssep,Tlong,
                                operator&=,
                                "Compute bitwise AND of image%s by image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_andeq,
                                "Compute bitwise AND of image%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute sequential bitwise AND of image%s");

        // Arctangent (two arguments).
        if (!std::strcmp("atan2",command)) {
          gmic_substitute_args(true);
          sep = 0;
          if (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                          gmic_use_indices,&sep,&end)==2 && sep==']' &&
              (ind=selection2cimg(indices,images.size(),images_names,"atan2")).height()==1) {
            print(images,0,"Compute pointwise oriented arctangent of image%s, "
                  "with x-argument [%u].",
                  gmic_selection.data(),
                  *ind);
            const CImg<T> img0 = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply(atan2(img0),true);
          } else arg_error("atan2");
          is_change = true;
          ++position;
          continue;
        }

        // Arccosine.
        gmic_simple_command("acos",acos,"Compute pointwise arccosine of image%s.");

        // Arcsine.
        gmic_simple_command("asin",asin,"Compute pointwise arcsine of image%s.");

        // Arctangent.
        gmic_simple_command("atan",atan,"Compute pointwise arctangent of image%s.");

        // Hyperbolic arccosine.
        gmic_simple_command("acosh",acosh,"Compute pointwise hyperbolic arccosine of image%s.");

        // Hyperbolic arcsine.
        gmic_simple_command("asinh",asinh,"Compute pointwise hyperbolic arcsine of image%s.");

        // Hyperbolic arctangent.
        gmic_simple_command("atanh",atanh,"Compute pointwise hyperbolic arctangent of image%s.");

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'b...'
        //-----------------------------
      gmic_commands_b :
        if (!is_get && item[1]=='r' && item[2]=='e' && item[3]=='a' && item[4]=='k' && !item[5]) // Redirect 'break'
          goto gmic_commands_others;

        // Blur.
        if (!std::strcmp("blur",command)) {
          gmic_substitute_args(false);
          unsigned int is_gaussian = 1;
          float sigma = -1;
          sep = *argx = 0;
          boundary = 1;

          const char *p_argument = argument;
          if (cimg_sscanf(argument,"%255[xyzc]%c",gmic_use_argx,&sep)==2 && sep==',') {
            p_argument+=1 + std::strlen(argx);
          } else sep = *argx = 0;

          if ((cimg_sscanf(p_argument,"%f%c",
                           &sigma,&end)==1 ||
               (cimg_sscanf(p_argument,"%f%c%c",
                            &sigma,&sep,&end)==2 && sep=='%') ||
               cimg_sscanf(p_argument,"%f,%u%c",
                           &sigma,&boundary,&end)==2 ||
               (cimg_sscanf(p_argument,"%f%c,%u%c",
                            &sigma,&sep,&boundary,&end)==3 && sep=='%') ||
               cimg_sscanf(p_argument,"%f,%u,%u%c",
                           &sigma,&boundary,&is_gaussian,&end)==3 ||
               (cimg_sscanf(p_argument,"%f%c,%u,%u%c",
                            &sigma,&sep,&boundary,&is_gaussian,&end)==4 && sep=='%')) &&
              sigma>=0 && boundary<=3 && is_gaussian<=1) {
            print(images,0,"Blur image%s%s%s%s with standard deviation %g%s, %s boundary conditions "
                  "and %s kernel.",
                  gmic_selection.data(),
                  *argx?" along the '":"",
                  *argx?argx:"",
                  *argx?std::strlen(argx)==1?"'-axis,":"'-axes,":"",
                  sigma,sep=='%'?"%":"",
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror",
                  is_gaussian?"gaussian":"quasi-gaussian");
            if (sep=='%') sigma = -sigma;
            if (*argx) {
              g_img.assign(4,1,1,1,(T)0);
              for (const char *s = argx; *s; ++s) g_img[*s>='x'?*s - 'x':3]+=sigma;
              cimg_forY(selection,l) gmic_apply(gmic_blur(g_img[0],g_img[1],g_img[2],g_img[3],
                                                          boundary,(bool)is_gaussian),true);
              g_img.assign();
            } else cimg_forY(selection,l) gmic_apply(blur(sigma,boundary,(bool)is_gaussian),true);
          } else arg_error("blur");
          is_change = true;
          ++position;
          continue;
        }

        // Box filter.
        if (!std::strcmp("boxfilter",command)) {
          unsigned int order = 0;
          gmic_substitute_args(false);
          float sigma = -1;
          sep = *argx = 0;
          boundary = 1;
          value = 1;
          const char *p_argument = argument;
          if (cimg_sscanf(argument,"%255[xyzc]%c",gmic_use_argx,&sep)==2 && sep==',') {
            p_argument+=1 + std::strlen(argx);
          } else sep = *argx = 0;
          if ((cimg_sscanf(p_argument,"%f%c",
                           &sigma,&end)==1 ||
               (cimg_sscanf(p_argument,"%f%c%c",
                            &sigma,&sep,&end)==2 && sep=='%') ||
               cimg_sscanf(p_argument,"%f,%u%c",
                           &sigma,&order,&end)==2 ||
               (cimg_sscanf(p_argument,"%f%c,%u%c",
                            &sigma,&sep,&order,&end)==3 && sep=='%') ||
               cimg_sscanf(p_argument,"%f,%u,%u%c",
                           &sigma,&order,&boundary,&end)==3 ||
               (cimg_sscanf(p_argument,"%f%c,%u,%u%c",
                            &sigma,&sep,&order,&boundary,&end)==4 && sep=='%') ||
               cimg_sscanf(p_argument,"%f,%u,%u,%lf%c",
                           &sigma,&order,&boundary,&value,&end)==4 ||
               (cimg_sscanf(p_argument,"%f%c,%u,%u,%lf%c",
                            &sigma,&sep,&order,&boundary,&value,&end)==5 && sep=='%')) &&
              sigma>=0 && boundary<=3 && order<=2 && value>=0) {
            const unsigned int nb_iter = (unsigned int)cimg::round(value);
            print(images,0,"Blur image%s%s%s%s with normalized box filter of size %g%s, order %u and "
                  "%s boundary conditions (%u iteration%s).",
                  gmic_selection.data(),
                  *argx?" along the '":"",
                  *argx?argx:"",
                  *argx?std::strlen(argx)==1?"'-axis,":"'-axes,":"",
                  sigma,sep=='%'?"%":"",
                  order,
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror",
                  nb_iter,nb_iter!=1?"s":"");
            if (sep=='%') sigma = -sigma;
            if (*argx) {
              g_img.assign(4,1,1,1,(T)0);
              for (const char *s = argx; *s; ++s) g_img[*s>='x'?*s - 'x':3]+=sigma;
              cimg_forY(selection,l) gmic_apply(gmic_blur_box(g_img[0],g_img[1],g_img[2],g_img[3],
                                                              order,boundary,nb_iter),true);
              g_img.assign();
            } else cimg_forY(selection,l) gmic_apply(gmic_blur_box(sigma,order,boundary,nb_iter),true);
          } else arg_error("boxfilter");
          is_change = true;
          ++position;
          continue;
        }

        // Bitwise right shift.
        gmic_arithmetic_command("bsr",
                                operator>>=,
                                "Compute bitwise right shift of image%s by %g%s",
                                gmic_selection.data(),value,ssep,Tlong,
                                operator>>=,
                                "Compute bitwise right shift of image%s by image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_brseq,
                                "Compute bitwise right shift of image%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute sequential bitwise right shift of image%s");

        // Bitwise left shift.
        gmic_arithmetic_command("bsl",
                                operator<<=,
                                "Compute bitwise left shift of image%s by %g%s",
                                gmic_selection.data(),value,ssep,Tlong,
                                operator<<=,
                                "Compute bitwise left shift of image%s by image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_blseq,
                                "Compute bitwise left shift of image%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute sequential bitwise left shift of image%s");

        // Bilateral filter.
        if (!std::strcmp("bilateral",command)) {
          gmic_substitute_args(true);
          float sigma_s = 0, sigma_r = 0, sampling_s = 0, sampling_r = 0;
          sep0 = sep1 = *argx = *argy = 0;
          if ((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           gmic_use_indices,gmic_use_argx,gmic_use_argy,&end)==3 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],%f,%f%c",
                           indices,argx,argy,&sampling_s,&sampling_r,&end)==5) &&
              (cimg_sscanf(argx,"%f%c",&sigma_s,&end)==1 ||
               (cimg_sscanf(argx,"%f%c%c",&sigma_s,&sep0,&end)==2 && sep0=='%')) &&
              (cimg_sscanf(argy,"%f%c",&sigma_r,&end)==1 ||
               (cimg_sscanf(argy,"%f%c%c",&sigma_r,&sep1,&end)==2 && sep1=='%')) &&
              (ind=selection2cimg(indices,images.size(),images_names,"bilateral")).height()==1 &&
              sigma_s>=0 && sigma_r>=0 && sampling_s>=0 && sampling_r>=0) {
            print(images,0,"Apply joint bilateral filter on image%s, with guide image [%u], "
                  " standard deviations (%g%s,%g%s) and sampling (%g,%g).",
                  gmic_selection.data(),
                  *ind,
                  sigma_s,sep0=='%'?"%":"",
                  sigma_r,sep1=='%'?"%":"",
                  sampling_s,sampling_r);
            const CImg<T> guide = gmic_image_arg(*ind);
            if (sep0=='%') sigma_s = -sigma_s;
            if (sep1=='%') sigma_r = -sigma_r;
            cimg_forY(selection,l) gmic_apply(blur_bilateral(guide,sigma_s,sigma_r,sampling_s,sampling_r),true);
          } else if ((cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                                  gmic_use_argx,gmic_use_argy,&end)==2 ||
                      cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%f,%f%c",
                                  argx,argy,&sampling_s,&sampling_r,&end)==4) &&
                     (cimg_sscanf(argx,"%f%c",&sigma_s,&end)==1 ||
                      (cimg_sscanf(argx,"%f%c%c",&sigma_s,&sep0,&end)==2 && sep0=='%')) &&
                     (cimg_sscanf(argy,"%f%c",&sigma_r,&end)==1 ||
                      (cimg_sscanf(argy,"%f%c%c",&sigma_r,&sep1,&end)==2 && sep1=='%')) &&
                     sigma_s>=0 && sigma_r>=0 && sampling_s>=0 && sampling_r>=0) {
            print(images,0,"Apply bilateral filter on image%s, with standard deviations (%g%s,%g%s) and "
                  "sampling (%g,%g).",
                  gmic_selection.data(),
                  sigma_s,sep0=='%'?"%":"",
                  sigma_r,sep1=='%'?"%":"",
                  sampling_s,sampling_r);
            if (sep0=='%') sigma_s = -sigma_s;
            if (sep1=='%') sigma_r = -sigma_r;
            cimg_forY(selection,l)
              gmic_apply(blur_bilateral(images[selection[l]],sigma_s,sigma_r,sampling_s,sampling_r),true);
          } else arg_error("bilateral");
          is_change = true;
          ++position;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'c...'
        //-----------------------------
      gmic_commands_c :

        // Check expression or filename.
        if (is_command_check) {
          gmic_substitute_args(false);
          is_cond = check_cond(argument,images,"check");
          if (is_very_verbose)
            print(images,0,"Check condition '%s' -> %s.",gmic_argument_text_printed(),is_cond?"true":"false");
          if (!is_cond) {
            if (is_first_item && callstack.size()>1 && callstack.back()[0]!='*')
              error(true,images,0,callstack.back(),"Command '%s': Invalid argument '%s'.",
                    callstack.back().data(),_gmic_argument_text(parent_arguments,gmic_use_argument_text,true));
            else error(true,images,0,0,
                       "Command 'check': Expression '%s' is false.",
                       gmic_argument_text());
          }
          ++position;
          continue;
        }

        // Crop.
        if (!std::strcmp("crop",command)) {
          gmic_substitute_args(false);
          name.assign(64,8);
          char
            *const st0 = name.data(0,0), *const st1 = name.data(0,1), *const st2 = name.data(0,2),
            *const st3 = name.data(0,3), *const st4 = name.data(0,4), *const st5 = name.data(0,5),
            *const st6 = name.data(0,6), *const st7 = name.data(0,7);
          char sep2 = sep0 = sep1 = 0, sep3 = 0, sep4 = 0, sep5 = 0, sep6 = 0, sep7 = 0;
          double a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0, a7 = 0;
          *st0 = *st1 = *st2 = *st3 = *st4 = *st5 = *st6 = *st7 = 0;
          boundary = 0;
          if ((boundary=0,cimg_sscanf(argument,"%63[0-9.eE%+-],%63[0-9.eE%+-]%c",
                                      st0,
                                      st1,&end)==2 ||
               cimg_sscanf(argument,"%63[0-9.eE%+-],%63[0-9.eE%+-],%u%c",
                           st0,
                           st1,&boundary,&end)==3) &&
              (cimg_sscanf(st0,"%lf%c",&a0,&end)==1 ||
               (cimg_sscanf(st0,"%lf%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
              (cimg_sscanf(st1,"%lf%c",&a1,&end)==1 ||
               (cimg_sscanf(st1,"%lf%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
              boundary<=3) {
            print(images,0,"Crop image%s with coordinates (%g%s) - (%g%s) and "
                  "%s boundary conditions.",
                  gmic_selection.data(),
                  a0,sep0=='%'?"%":"",
                  a1,sep1=='%'?"%":"",
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*(img.width() - 1)/100:a0),
                x1 = (int)cimg::round(sep1=='%'?a1*(img.width() - 1)/100:a1);
              if (img) { gmic_apply(crop(x0,x1,boundary),false); }
              else error(true,images,0,0,
                         "Command 'crop': Cannot crop empty image [%d] (with arguments '%s').",
                         selection[l],gmic_argument_text());
            }
          } else if ((boundary=0,cimg_sscanf(argument,
                                             "%63[0-9.eE%+-],%63[0-9.eE%+-],"
                                             "%63[0-9.eE%+-],%63[0-9.eE%+-]%c",
                                             st0,st1,st2,st3,&end)==4 ||
                      cimg_sscanf(argument,
                                  "%63[0-9.eE%+-],%63[0-9.eE%+-],"
                                  "%63[0-9.eE%+-],%63[0-9.eE%+-],%u%c",
                                  st0,st1,st2,st3,&boundary,&end)==5) &&
                     (cimg_sscanf(st0,"%lf%c",&a0,&end)==1 ||
                      (cimg_sscanf(st0,"%lf%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (cimg_sscanf(st1,"%lf%c",&a1,&end)==1 ||
                      (cimg_sscanf(st1,"%lf%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
                     (cimg_sscanf(st2,"%lf%c",&a2,&end)==1 ||
                      (cimg_sscanf(st2,"%lf%c%c",&a2,&sep2,&end)==2 && sep2=='%')) &&
                     (cimg_sscanf(st3,"%lf%c",&a3,&end)==1 ||
                      (cimg_sscanf(st3,"%lf%c%c",&a3,&sep3,&end)==2 && sep3=='%')) &&
                     boundary<=3) {
            print(images,0,
                  "Crop image%s with coordinates (%g%s,%g%s) - (%g%s,%g%s) and "
                  "%s boundary conditions.",
                  gmic_selection.data(),
                  a0,sep0=='%'?"%":"",
                  a1,sep1=='%'?"%":"",
                  a2,sep2=='%'?"%":"",
                  a3,sep3=='%'?"%":"",
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*(img.width() - 1)/100:a0),
                y0 = (int)cimg::round(sep1=='%'?a1*(img.height() - 1)/100:a1),
                x1 = (int)cimg::round(sep2=='%'?a2*(img.width() - 1)/100:a2),
                y1 = (int)cimg::round(sep3=='%'?a3*(img.height() - 1)/100:a3);
              if (img) { gmic_apply(crop(x0,y0,x1,y1,boundary),false); }
              else error(true,images,0,0,
                         "Command 'crop': Cannot crop empty image [%d] (with arguments '%s').",
                         selection[l],gmic_argument_text());
            }
          } else if ((boundary=0,cimg_sscanf(argument,
                                             "%63[0-9.eE%+-],%63[0-9.eE%+-],%63[0-9.eE%+-],"
                                             "%63[0-9.eE%+-],%63[0-9.eE%+-],%63[0-9.eE%+-]%c",
                                             st0,st1,st2,st3,st4,st5,&end)==6 ||
                      cimg_sscanf(argument,"%63[0-9.eE%+-],%63[0-9.eE%+-],%63[0-9.eE%+-],"
                                  "%63[0-9.eE%+-],%63[0-9.eE%+-],%63[0-9.eE%+-],%u%c",
                                  st0,st1,st2,st3,st4,st5,&boundary,&end)==7) &&
                     (cimg_sscanf(st0,"%lf%c",&a0,&end)==1 ||
                      (cimg_sscanf(st0,"%lf%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (cimg_sscanf(st1,"%lf%c",&a1,&end)==1 ||
                      (cimg_sscanf(st1,"%lf%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
                     (cimg_sscanf(st2,"%lf%c",&a2,&end)==1 ||
                      (cimg_sscanf(st2,"%lf%c%c",&a2,&sep2,&end)==2 && sep2=='%')) &&
                     (cimg_sscanf(st3,"%lf%c",&a3,&end)==1 ||
                      (cimg_sscanf(st3,"%lf%c%c",&a3,&sep3,&end)==2 && sep3=='%')) &&
                     (cimg_sscanf(st4,"%lf%c",&a4,&end)==1 ||
                      (cimg_sscanf(st4,"%lf%c%c",&a4,&sep4,&end)==2 && sep4=='%')) &&
                     (cimg_sscanf(st5,"%lf%c",&a5,&end)==1 ||
                      (cimg_sscanf(st5,"%lf%c%c",&a5,&sep5,&end)==2 && sep5=='%')) &&
                     boundary<=3) {
            print(images,0,"Crop image%s with coordinates (%g%s,%g%s,%g%s) - (%g%s,%g%s,%g%s) "
                  "and %s boundary conditions.",
                  gmic_selection.data(),
                  a0,sep0=='%'?"%":"",
                  a1,sep1=='%'?"%":"",
                  a2,sep2=='%'?"%":"",
                  a3,sep3=='%'?"%":"",
                  a4,sep4=='%'?"%":"",
                  a5,sep5=='%'?"%":"",
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*(img.width() - 1)/100:a0),
                y0 = (int)cimg::round(sep1=='%'?a1*(img.height() - 1)/100:a1),
                z0 = (int)cimg::round(sep2=='%'?a2*(img.depth() - 1)/100:a2),
                x1 = (int)cimg::round(sep3=='%'?a3*(img.width() - 1)/100:a3),
                y1 = (int)cimg::round(sep4=='%'?a4*(img.height() - 1)/100:a4),
                z1 = (int)cimg::round(sep5=='%'?a5*(img.depth() - 1)/100:a5);
              if (img) { gmic_apply(crop(x0,y0,z0,x1,y1,z1,boundary),false); }
              else error(true,images,0,0,
                         "Command 'crop': Cannot crop empty image [%d] (with arguments '%s').",
                         selection[l],gmic_argument_text());
            }
          } else if ((boundary=0,cimg_sscanf(argument,
                                             "%63[0-9.eE%+-],%63[0-9.eE%+-],%63[0-9.eE%+-],"
                                             "%63[0-9.eE%+-],%63[0-9.eE%+-],%63[0-9.eE%+-],"
                                             "%63[0-9.eE%+-],%63[0-9.eE%+-]%c",
                                             st0,st1,st2,st3,st4,st5,st6,st7,&end)==8 ||
                      cimg_sscanf(argument,"%63[0-9.eE%+-],%63[0-9.eE%+-],%63[0-9.eE%+-],"
                                  "%63[0-9.eE%+-],%63[0-9.eE%+-],%63[0-9.eE%+-],"
                                  "%63[0-9.eE%+-],%63[0-9.eE%+-],%u%c",
                                  st0,st1,st2,st3,st4,st5,st6,st7,&boundary,&end)==9) &&
                     (cimg_sscanf(st0,"%lf%c",&a0,&end)==1 ||
                      (cimg_sscanf(st0,"%lf%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (cimg_sscanf(st1,"%lf%c",&a1,&end)==1 ||
                      (cimg_sscanf(st1,"%lf%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
                     (cimg_sscanf(st2,"%lf%c",&a2,&end)==1 ||
                      (cimg_sscanf(st2,"%lf%c%c",&a2,&sep2,&end)==2 && sep2=='%')) &&
                     (cimg_sscanf(st3,"%lf%c",&a3,&end)==1 ||
                      (cimg_sscanf(st3,"%lf%c%c",&a3,&sep3,&end)==2 && sep3=='%')) &&
                     (cimg_sscanf(st4,"%lf%c",&a4,&end)==1 ||
                      (cimg_sscanf(st4,"%lf%c%c",&a4,&sep4,&end)==2 && sep4=='%')) &&
                     (cimg_sscanf(st5,"%lf%c",&a5,&end)==1 ||
                      (cimg_sscanf(st5,"%lf%c%c",&a5,&sep5,&end)==2 && sep5=='%')) &&
                     (cimg_sscanf(st6,"%lf%c",&a6,&end)==1 ||
                      (cimg_sscanf(st6,"%lf%c%c",&a6,&sep6,&end)==2 && sep6=='%')) &&
                     (cimg_sscanf(st7,"%lf%c",&a7,&end)==1 ||
                      (cimg_sscanf(st7,"%lf%c%c",&a7,&sep7,&end)==2 && sep7=='%')) &&
                     boundary<=3) {
            print(images,0,
                  "Crop image%s with coordinates (%g%s,%g%s,%g%s,%g%s) - (%g%s,%g%s,%g%s,%g%s) "
                  "and %s boundary conditions.",
                  gmic_selection.data(),
                  a0,sep0=='%'?"%":"",
                  a1,sep1=='%'?"%":"",
                  a2,sep2=='%'?"%":"",
                  a3,sep3=='%'?"%":"",
                  a4,sep4=='%'?"%":"",
                  a5,sep5=='%'?"%":"",
                  a6,sep6=='%'?"%":"",
                  a7,sep7=='%'?"%":"",
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                x0 = (int)cimg::round(sep0=='%'?a0*(img.width() - 1)/100:a0),
                y0 = (int)cimg::round(sep1=='%'?a1*(img.height() - 1)/100:a1),
                z0 = (int)cimg::round(sep2=='%'?a2*(img.depth() - 1)/100:a2),
                v0 = (int)cimg::round(sep3=='%'?a3*(img.spectrum() - 1)/100:a3),
                x1 = (int)cimg::round(sep4=='%'?a4*(img.width() - 1)/100:a4),
                y1 = (int)cimg::round(sep5=='%'?a5*(img.height() - 1)/100:a5),
                z1 = (int)cimg::round(sep6=='%'?a6*(img.depth() - 1)/100:a6),
                v1 = (int)cimg::round(sep7=='%'?a7*(img.spectrum() - 1)/100:a7);
              if (img) { gmic_apply(crop(x0,y0,z0,v0,x1,y1,z1,v1,boundary),false); }
              else error(true,images,0,0,
                         "Command 'crop': Cannot crop empty image [%d] (with arguments '%s').",
                         selection[l],gmic_argument_text());
            }
          } else arg_error("crop");
          is_change = true;
          ++position;
          continue;
        }

        // Cut.
        if (!std::strcmp("cut",command)) {
          gmic_substitute_args(true);
          ind0.assign(); ind1.assign();
          sep0 = sep1 = *argx = *argy = *indices = 0;
          value0 = value1 = 0;
          if (cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-]%c",
                          gmic_use_argx,gmic_use_argy,&end)==2 &&
              ((cimg_sscanf(argx,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep0,&end)==2 &&
                sep0==']' &&
                (ind0=selection2cimg(indices,images.size(),images_names,"cut")).height()==1) ||
               (cimg_sscanf(argx,"%lf%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
               cimg_sscanf(argx,"%lf%c",&value0,&end)==1) &&
              ((cimg_sscanf(argy,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_formula,&sep1,&end)==2 &&
                sep1==']' &&
                (ind1=selection2cimg(formula,images.size(),images_names,"cut")).height()==1) ||
               (cimg_sscanf(argy,"%lf%c%c",&value1,&sep1,&end)==2 && sep1=='%') ||
               cimg_sscanf(argy,"%lf%c",&value1,&end)==1)) {
            if (ind0) { value0 = images[*ind0].min(); sep0 = 0; }
            if (ind1) { value1 = images[*ind1].max(); sep1 = 0; }
            print(images,0,"Cut image%s in range [%g%s,%g%s].",
                  gmic_selection.data(),
                  value0,sep0=='%'?"%":"",
                  value1,sep1=='%'?"%":"");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              nvalue0 = value0; nvalue1 = value1;
              vmin = vmax = 0;
              if (sep0=='%' || sep1=='%') {
                if (img) vmax = (double)img.max_min(vmin);
                if (sep0=='%') nvalue0 = vmin + (vmax - vmin)*value0/100;
                if (sep1=='%') nvalue1 = vmin + (vmax - vmin)*value1/100;
              }
              gmic_apply(cut((T)nvalue0,(T)nvalue1),true);
            }
          } else if (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep0,&end)==2 &&
                     sep0==']' &&
                     (ind0=selection2cimg(indices,images.size(),images_names,"cut")).height()==1) {
            if (images[*ind0]) value1 = (double)images[*ind0].max_min(value0);
            print(images,0,"Cut image%s in range [%g,%g].",
                  gmic_selection.data(),
                  value0,
                  value1);
            cimg_forY(selection,l) gmic_apply(cut((T)value0,(T)value1),true);
          } else arg_error("cut");
          is_change = true;
          ++position;
          continue;
        }

        // Import custom commands.
        if (!is_get && !std::strcmp("command",item)) {
          gmic_substitute_args(false);
          name.assign(argument,(unsigned int)std::strlen(argument) + 1);
          const char *arg_command_text = gmic_argument_text_printed();
          unsigned int offset_argument_text = 0, count_new = 0, count_replaced = 0;
          char *arg_command = name;
          strreplace_fw(arg_command);

          bool add_debug_info = is_debug;
          const bool is_debug_arg = (*arg_command=='0' || *arg_command=='1') && arg_command[1]==',';
          if (is_debug_arg) {
            add_debug_info = (*arg_command=='1');
            arg_command+=2; arg_command_text+=2; offset_argument_text = 2;
          }

          std::FILE *file = 0;
          const char *const p_colon = std::strchr(arg_command,':');
#if cimg_OS!=2
          if (!p_colon || p_colon - arg_command>256)
            file = cimg::std_fopen(arg_command,"rb");
#else
          if (!p_colon || p_colon - arg_command==1 || p_colon - arg_command>256)
            file = cimg::std_fopen(arg_command,"rb"); // Allow 'C:\\filename'
#endif
          if (file) {
            if (!is_debug_arg) add_debug_info = true;
            print(images,0,"Import commands from file '%s'%s",
                  arg_command_text,
                  add_debug_info?", with debug info":"");
            add_commands(file,arg_command,add_debug_info,&count_new,&count_replaced);
            cimg::fclose(file);

          } else if (!cimg::strncasecmp(arg_command,"http://",7) ||
                     !cimg::strncasecmp(arg_command,"https://",8)) { // Try to read from network
            if (!is_debug_arg) add_debug_info = true;
            print(images,0,"Import commands from URL '%s'%s",
                  arg_command_text,
                  add_debug_info?", with debug info":"");
            try {
              file = cimg::std_fopen(cimg::load_network(arg_command,gmic_use_argx,network_timeout,true,0,"gmic"),"r");
            } catch (...) {
              file = 0;
            }
            if (file) {
              status.move_to(o_status); // Save status because 'add_commands()' can change it, with 'error()'
              const int o_verbosity = verbosity;
              const bool o_is_debug = is_debug;
              verbosity = 0;
              is_debug = false;
              try {
                add_commands(file,arg_command,add_debug_info,&count_new,&count_replaced);
                cimg::fclose(file);
              } catch (...) {
                cimg::fclose(file);
                file = 0;
              }
              is_debug = o_is_debug;
              verbosity = o_verbosity;
              o_status.move_to(status);
            }
            if (!file)
              error(true,images,0,0,
                    "Command 'command': Unable to load custom command file '%s' "
                    "from network.",
                    gmic_argument_text() + offset_argument_text);
            std::remove(argx);
          } else { // Import commands from a string
            if (!is_debug_arg) add_debug_info = false;
            print(images,0,"Import custom commands from expression '%s'%s",
                  arg_command_text,
                  add_debug_info?", with debug info":"");
            cimg::strunescape(arg_command);
            add_commands(arg_command,0,add_debug_info,&count_new,&count_replaced);
          }
          if (is_verbose) {
            unsigned int count_total = 0;
            for (unsigned int l = 0; l<gmic_comslots; ++l) count_total+=commands[l].size();
            cimg::mutex(29);
            if (count_new && count_replaced)
              std::fprintf(cimg::output()," (%u new, %u replaced, total: %u).",
                           count_new,count_replaced,count_total);
            else if (count_new)
              std::fprintf(cimg::output()," (%u new, total: %u).",
                           count_new,count_total);
            else
              std::fprintf(cimg::output()," (%u replaced, total: %u).",
                           count_replaced,count_total);
            std::fflush(cimg::output());
            cimg::mutex(29,0);
          }
          ++position;
          continue;
        }

        // Check validity of 3D object.
        if (!is_get && !std::strcmp("check3d",command)) {
          gmic_substitute_args(false);
          bool is_full_check = false;
          if ((*argument=='0' || *argument=='1') && !argument[1]) {
            is_full_check = (*argument=='1');
            ++position;
          } else is_full_check = false;
          if (is_very_verbose)
            print(images,0,"Check validity of 3D object%s (%s check)",
                  gmic_selection.data(),
                  is_full_check?"full":"fast");
          cimg_forY(selection,l) {
            uind = selection[l];
            CImg<T>& img = gmic_check(images[uind]);
            if (!img.is_CImg3d(is_full_check,&(*gmic_use_message=0))) {
              if (is_very_verbose) {
                cimg::mutex(29);
                std::fprintf(cimg::output()," -> invalid.");
                std::fflush(cimg::output());
                cimg::mutex(29,0);
              }
              error(true,images,0,0,
                    "Command 'check3d': Invalid 3D object [%d], in image%s (%s).",
                    uind,gmic_selection_err.data(),message);
            }
          }
          if (is_very_verbose) {
            cimg::mutex(29);
            std::fprintf(cimg::output()," -> valid.");
            std::fflush(cimg::output());
            cimg::mutex(29,0);
          }
          continue;
        }

        // Cosine.
        gmic_simple_command("cos",cos,"Compute pointwise cosine of image%s.");

        // Convolve & Correlate.
        if (!std::strcmp("convolve",command) || !std::strcmp("correlate",command)) {
          gmic_substitute_args(true);
          unsigned int is_normalized = 0, channel_mode = 1, interpolation_type = 0;
          int
            xstart = 0, ystart = 0, zstart = 0,
            xend = (int)(~0U>>1), yend = (int)(~0U>>1), zend = (int)(~0U>>1),
            xcenter = (int)(~0U>>1), ycenter = (int)(~0U>>1), zcenter = (int)(~0U>>1);
          float
            xstride = 1, ystride = 1, zstride = 1,
            xdilation = 1, ydilation = 1 , zdilation = 1;
          is_cond = command[2]=='n'; // is_convolve?
          boundary = 1;
          sep = 0;

          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                            gmic_use_indices,&sep,&end)==2 && sep==']') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u%c",
                           indices,&boundary,&end)==2 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u%c",
                           indices,&boundary,&is_normalized,&end)==3 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%u%c",
                           indices,&boundary,&is_normalized,&channel_mode,&end)==4 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%u,%d,%d,%d%c",
                           indices,&boundary,&is_normalized,&channel_mode,&xcenter,&ycenter,&zcenter,&end)==7 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%u,%d,%d,%d,%d,%d,%d,%d,%d,%d%c",
                           indices,&boundary,&is_normalized,&channel_mode,&xcenter,&ycenter,&zcenter,
                           &xstart,&ystart,&zstart,&xend,&yend,&zend,&end)==13 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%u,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f%c",
                           indices,&boundary,&is_normalized,&channel_mode,&xcenter,&ycenter,&zcenter,
                           &xstart,&ystart,&zstart,&xend,&yend,&zend,&xstride,&ystride,&zstride,&end)==16 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%u,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f%c",
                           indices,&boundary,&is_normalized,&channel_mode,&xcenter,&ycenter,&zcenter,
                           &xstart,&ystart,&zstart,&xend,&yend,&zend,&xstride,&ystride,&zstride,
                           &xdilation,&ydilation,&zdilation,&end)==19 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%u,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%u%c",
                           indices,&boundary,&is_normalized,&channel_mode,&xcenter,&ycenter,&zcenter,
                           &xstart,&ystart,&zstart,&xend,&yend,&zend,&xstride,&ystride,&zstride,
                           &xdilation,&ydilation,&zdilation,&interpolation_type,&end)==20) &&
              (ind=selection2cimg(indices,images.size(),images_names,"correlate")).height()==1 &&
              boundary<=3 && channel_mode<=3 && interpolation_type<=1) {

            *argx = *argy = *argz = *argc = 0;
            if (is_verbose) {
              if (xcenter!=(int)(~0U>>1) || ycenter!=(int)(~0U>>1) || zcenter!=(int)(~0U>>1)) {
                gmic_use_argx;
                cimg_snprintf(argx,_argx.width(),", kernel center (%d,%d,%d)",
                              (int)xcenter,(int)ycenter,(int)zcenter);
              }
              if (xstart!=0 || ystart!=0 || zstart!=0 ||
                  xend!=(int)(~0U>>1) || yend!=(int)(~0U>>1) || zend!=(int)(~0U>>1)) {
                gmic_use_argy;
                cimg_snprintf(argy,_argy.width(),", crop coordinates (%d,%d,%d) - (%d,%d,%d)",
                              xstart,ystart,zstart,xend,yend,zend);
              }
              if (xstride!=1 || ystride!=1 || zstride!=1) {
                gmic_use_argz;
                cimg_snprintf(argz,_argz.width(),", strides (%g,%g,%g)",
                              xstride,ystride,zstride);
              }
              if (xdilation!=1 || ydilation!=1 || zdilation!=1) {
                gmic_use_argc;
                cimg_snprintf(argc,_argc.width(),", dilations (%g,%g,%g)",
                              xdilation,ydilation,zdilation);
              }
            }

            print(images,0,
                  "%s image%s with kernel [%u], %s boundary conditions, "
                  "with%s normalization, channel mode '%s'%s%s%s%s and %s interpolation.",
                  is_cond?"Convolve":"Correlate",
                  gmic_selection.data(),
                  *ind,
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror",
                  is_normalized?"":"out",
                  channel_mode==0?"all":
                  channel_mode==1?"one for one":
                  channel_mode==2?"partial sum":"full sum",
                  *argx?argx:"",*argy?argy:"",*argz?argz:"",*argc?argc:"",
                  interpolation_type?"linear":"nearest-neighbor");
            const CImg<T> kernel = gmic_image_arg(*ind);
            if (is_cond) {
              cimg_forY(selection,l) gmic_apply(convolve(kernel,boundary,(bool)is_normalized,channel_mode,
                                                         xcenter,ycenter,zcenter,xstart,ystart,zstart,xend,yend,zend,
                                                         xstride,ystride,zstride,xdilation,ydilation,zdilation,
                                                         interpolation_type),false);
            } else {
              cimg_forY(selection,l) gmic_apply(correlate(kernel,boundary,(bool)is_normalized,channel_mode,
                                                          xcenter,ycenter,zcenter,xstart,ystart,zstart,xend,yend,zend,
                                                          xstride,ystride,zstride,xdilation,ydilation,zdilation,
                                                          interpolation_type),false);
            }
          } else arg_error(is_cond?"convolve":"correlate");
          is_change = true;
          ++position;
          continue;
        }

        // Cumulate.
        if (!std::strcmp("cumulate",command)) {
          gmic_substitute_args(false);
          bool is_axes_argument = true;
          for (const char *s = argument; *s; ++s) {
            const char _s = *s;
            if (_s!='x' && _s!='y' && _s!='z' && _s!='c') { is_axes_argument = false; break; }
          }
          if (*argument && is_axes_argument) {
            print(images,0,"Cumulate values of image%s along the '%s'-ax%cs.",
                  gmic_selection.data(),
                  gmic_argument_text_printed(),
                  std::strlen(argument)>1?'e':'i');
            cimg_forY(selection,l) gmic_apply(cumulate(argument),true);
            ++position;
          } else {
            print(images,0,"Cumulate values of image%s.",
                  gmic_selection.data());
            cimg_forY(selection,l) gmic_apply(cumulate(),true);
          }
          is_change = true;
          continue;
        }

        // Hyperbolic cosine.
        gmic_simple_command("cosh",cosh,"Compute pointwise hyperbolic cosine of image%s.");

        // Camera input.
        if (!std::strcmp("camera",item)) {
          gmic_substitute_args(false);
          double
            cam_index = 0, nb_frames = 1, skip_frames = 0,
            capture_width = 0, capture_height = 0;
          if ((cimg_sscanf(argument,"%lf%c",
                           &cam_index,&end)==1 ||
               cimg_sscanf(argument,"%lf,%lf%c",
                           &cam_index,&nb_frames,&end)==2 ||
               cimg_sscanf(argument,"%lf,%lf,%lf%c",
                           &cam_index,&nb_frames,&skip_frames,&end)==3 ||
               cimg_sscanf(argument,"%lf,%lf,%lf,%lf,%lf%c",
                           &cam_index,&nb_frames,&skip_frames,
                           &capture_width,&capture_height,&end)==5) &&
              cam_index>=0 && nb_frames>=0 && skip_frames>=0 &&
              ((!capture_width && !capture_height) || (capture_width>0 && capture_height>0)))
            ++position;
          else { cam_index = skip_frames = capture_width = capture_height = 0; nb_frames = 1; }
          cam_index = cimg::round(cam_index);
          nb_frames = cimg::round(nb_frames);
          skip_frames = cimg::round(skip_frames);
          capture_width = cimg::round(capture_width);
          capture_height = cimg::round(capture_height);
          if (!nb_frames) {
            print(images,0,"Release camera #%g.",cam_index);
            CImg<T>::get_load_camera((unsigned int)cam_index,0,0,0,true);
          } else {
            if (capture_width)
              print(images,0,"Insert %g image%s from camera #%g, with %g frames skipping "
                    "and resolution %gx%g.",
                    cam_index,nb_frames,nb_frames>1?"s":"",skip_frames,
                    capture_width,capture_height);
            else print(images,0,"Insert %g image%s from camera #%g, with %g frames skipping.",
                       cam_index,nb_frames,nb_frames>1?"s":"",skip_frames);
            gmic_use_title;
            cimg_snprintf(title,_title.width(),"[Camera #%g]",cam_index);
            CImg<char>::string(title).move_to(name);
            if (nb_frames>1) {
              cimg::mutex(29);
              std::fputc('\n',cimg::output());
              std::fflush(cimg::output());
              cimg::mutex(29,0);
            }
            for (unsigned int k = 0; k<(unsigned int)nb_frames; ++k) {
              if (nb_frames>1 && is_verbose) {
                cimg::mutex(29);
                std::fprintf(cimg::output(),"\r  > Image %u/%u        ",
                             k + 1,(unsigned int)nb_frames);
                std::fflush(cimg::output());
                cimg::mutex(29,0);
              }
              CImg<T>::get_load_camera((unsigned int)cam_index,
                                       (unsigned int)capture_width,(unsigned int)capture_height,
                                       (unsigned int)skip_frames,false).
                move_to(images);
              images_names.insert(name);
            }
          }
          is_change = true;
          continue;
        }

        // Show/hide mouse cursor.
        if (!is_get && !std::strcmp("cursor",command)) {
          gmic_substitute_args(false);
          if (!is_selection)
            CImg<unsigned int>::vector(0,1,2,3,4,5,6,7,8,9).move_to(selection);
          bool state = true;
          if ((*argument=='0' || *argument=='1') && !argument[1]) {
            state = (*argument=='1'); ++position;
          } else state = true;

          if (!is_display_available) {
            print(images,0,"%s mouse cursor for display window%s (skipped, no display %s).",
                  state?"Show":"Hide",
                  gmic_selection.data(),cimg_display?"available":"support");
          } else {
            if (state) cimg_forY(selection,l) gmic_display_window(selection[l]).show_mouse();
            else cimg_forY(selection,l) gmic_display_window(selection[l]).hide_mouse();
            print(images,0,"%s mouse cursor for display window%s.",
                  state?"Show":"Hide",
                  gmic_selection.data());
          }
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'd...'
        //-----------------------------
      gmic_commands_d :
        if (command[1]=='i' && command[2]=='v' && command[3]=='3' && command[4]=='d' && !command[5]) // Redirect 'div3d'
          goto gmic_commands_others;

        // Done.
        if (!is_get && !std::strcmp("done",item)) {
          const CImg<char> &s = callstack.back();
          if (s[0]!='*' || (s[1]!='b' && s[1]!='f' && s[1]!='l' && s[1]!='r'))
            error(true,images,0,0,
                  "Command 'done': Not associated to a 'for', 'foreach', 'local' or 'repeat' command "
                  "within the same scope.");

          if (s[1]=='b') callstack.remove(); // End a '{ .. }' block
          else if (s[1]=='f') {
            if (s[4]!='e') { // End a 'for...done' block
              unsigned int *const fd = fordones.data(0,nb_fordones - 1);
              position = fd[0];
              ++fd[1];
              next_debug_line = fd[2];
              next_debug_filename = debug_filename;
            } else { // End a 'foreach...done' block
              is_end_local = true;
              break;
            }
          } else if (s[1]=='l') { // End a 'local...done' block
            if (is_very_verbose) print(images,0,"End 'local...done' block.");
            is_end_local = true;
            break;
          } else if (s[1]=='r') { // End a 'repeat...done' block
            unsigned int *const rd = repeatdones.data(0,nb_repeatdones - 1);
            ++rd[1];
            if (rd[2]!=~0U) --rd[2];
            if (rd[2]) {
              position = rd[0] + 2;
              next_debug_line = rd[3];
              next_debug_filename = debug_filename;
              is_lbrace_command = true;
            } else {
              if (is_very_verbose) print(images,0,"End 'repeat...done' block.");
              --nb_repeatdones;
              callstack.remove();
            }
          }
          continue;
        }

        // Do...while.
        if (!is_get && !std::strcmp("do",item)) {
          if (is_debug_info && debug_line!=~0U) {
            gmic_use_argx;
            cimg_snprintf(argx,_argx.width(),"*do#%u",debug_line);
            CImg<char>::string(argx).move_to(callstack);
          } else CImg<char>::string("*do").move_to(callstack);
          if (is_very_verbose) print(images,0,"Start 'do...while' block.");
          if (nb_dowhiles>=dowhiles._height) dowhiles.resize(3,std::max(2*dowhiles._height,8U),1,1,0);
          unsigned int *const dw = dowhiles.data(0,nb_dowhiles++);
          dw[0] = position_item;
          dw[1] = 0;
          dw[2] = debug_line;
          continue;
        }

        // Discard value.
        if (!std::strcmp("discard",command)) {
          gmic_substitute_args(false);
          CImg<T> values;
          *argx = 0;

          if (cimg_sscanf(argument,"%255[xyzc]%c",gmic_use_argx,&end)==1) {

            // Discard neighboring duplicate values along axes.
            print(images,0,"Discard neighboring duplicate values along '%s'-ax%cs, in image%s.",
                  argx,
                  std::strlen(argx)>1?'e':'i',
                  gmic_selection.data());
            cimg_forY(selection,l) gmic_apply(gmic_discard(argx),false);
            ++position;

          } else if (cimg_sscanf(argument,"%255[xyzc]%c",gmic_use_argx,&end)==2 && end==',') {

            // Discard sequence of values along axes.
            unsigned int nb_values = 1, s_argx = (unsigned int)std::strlen(argx) + 1;
            for (const char *s = argument + s_argx; *s; ++s) if (*s==',') ++nb_values;
            try { values.assign(nb_values,1,1,1).fill_from_values(argument + s_argx,true); }
            catch (CImgException&) { arg_error("discard"); }
            print(images,0,"Discard sequence of values '%s' along '%s'-ax%cs, in image%s.",
                  gmic_argument_text_printed() + s_argx,
                  argx,
                  std::strlen(argx)>1?'e':'i',
                  gmic_selection.data());
            cimg_forY(selection,l) gmic_apply(gmic_discard(values,argx),false);
            ++position;

          } else { // Discard sequence of values or neighboring duplicate values
            unsigned int nb_values = *argument?1U:0U;
            for (const char *s = argument; *s; ++s) if (*s==',') ++nb_values;
            try { values.assign(nb_values,1,1,1).fill_from_values(argument,true); }
            catch (CImgException&) { values.assign(); }
            if (values) {
              print(images,0,"Discard sequence of values '%s' in image%s.",
                    gmic_argument_text_printed(),
                    gmic_selection.data());
              cimg_forY(selection,l) gmic_apply(discard(values),false);
              ++position;

            } else {
              print(images,0,"Discard neighboring duplicate values in image%s.",
                    gmic_selection.data());
              cimg_forY(selection,l) gmic_apply(discard(),false);
            }
          }
          is_change = true;
          continue;
        }

        // Enable debug mode (useful when 'debug' is invoked from a custom command).
        if (!is_get && !std::strcmp("debug",item)) {
          is_debug = true;
          continue;
        }

        // Divide.
        gmic_arithmetic_command("div",
                                operator/=,
                                "Divide image%s by %g%s",
                                gmic_selection.data(),value,ssep,Tfloat,
                                div,
                                "Divide image%s by image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_diveq,
                                "Divide image%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Divide image%s");

        // Distance function.
        if (!std::strcmp("distance",command)) {
          gmic_substitute_args(true);
          unsigned int algorithm = 0, off = 0;
          int metric = 2;
          sep0 = sep1 = *indices = 0;
          value = 0;
          if ((cimg_sscanf(argument,"%lf%c",
                           &value,&end)==1 ||
               (cimg_sscanf(argument,"%lf%c%c",
                            &value,&sep0,&end)==2 && sep0=='%') ||
               cimg_sscanf(argument,"%lf,%d%c",
                           &value,&metric,&end)==2 ||
               (cimg_sscanf(argument,"%lf%c,%d%c",
                            &value,&sep0,&metric,&end)==3 && sep0=='%')) &&
              metric>=0 && metric<=3) {
            print(images,0,"Compute distance map to isovalue %g%s in image%s, "
                  "with %s metric.",
                  value,sep0=='%'?"%":"",
                  gmic_selection.data(),
                  metric==0?"chebyshev":metric==1?"manhattan":metric==2?"euclidean":
                  "squared-euclidean");
            cimg_forY(selection,l) {
              CImg<T> &img = gmic_check(images[selection[l]]);
              nvalue = value;
              if (sep0=='%' && img) {
                vmax = (double)img.max_min(vmin);
                nvalue = vmin + value*(vmax - vmin)/100;
              }
              gmic_apply(distance((T)nvalue,metric),true);
            }
          } else if ((((cimg_sscanf(argument,"%lf,[%255[a-zA-Z0-9_.%+-]%c%c",
                                    &value,gmic_use_indices,&sep1,&end)==3 ||
                        (cimg_sscanf(argument,"%lf%c,[%255[a-zA-Z0-9_.%+-]%c%c",
                                     &value,&sep0,indices,&sep1,&end)==4 && sep0=='%')) &&
                       sep1==']') ||
                      ((cimg_sscanf(argument,"%lf,[%255[a-zA-Z0-9_.%+-]],%u%c",
                                    &value,indices,&algorithm,&end)==3 ||
                        (cimg_sscanf(argument,"%lf%c,[%255[a-zA-Z0-9_.%+-]],%u%c",
                                     &value,&sep0,indices,&algorithm,&end)==4 && sep0=='%')) &&
                       algorithm<=4)) &&
                     (ind=selection2cimg(indices,images.size(),images_names,"distance")).height()==1) {
            print(images,0,"Compute distance map%s to isovalue %g%s in image%s, "
                  "using %s algorithm, with metric [%u].",
                  selection.height()>1?(algorithm>=3?"s and return paths":"s"):
                  (algorithm>=3?" and return path":""),
                  value,sep0=='%'?"%":"",
                  gmic_selection.data(),
                  algorithm==0?"fast-marching":algorithm==1||algorithm==3?
                  "low-connectivity dijkstra":"high-connectivity dijkstra",
                  *ind);
            const CImg<T> custom_metric = gmic_image_arg(*ind);
            if (algorithm<3) cimg_forY(selection,l) {
                CImg<T> &img = gmic_check(images[selection[l]]);
                nvalue = value;
                if (sep0=='%' && img) {
                  vmax = (double)img.max_min(vmin);
                  nvalue = vmin + value*(vmax - vmin)/100;
                }
                if (!algorithm) { gmic_apply(distance_eikonal((T)nvalue,custom_metric),false); }
                else gmic_apply(distance_dijkstra((T)nvalue,custom_metric,algorithm==2),false);
              }
            else cimg_forY(selection,l) {
                uind = selection[l] + off;
                CImg<T>& img = gmic_check(images[uind]);
                nvalue = value;
                if (sep0=='%' && img) {
                  vmax = (double)img.max_min(vmin);
                  nvalue = vmin + value*(vmax - vmin)/100;
                }
                CImg<T> path(1),
                  dist = img.get_distance_dijkstra((T)nvalue,custom_metric,algorithm==4,path);
                dist.append(path,'c');
                if (is_get) {
                  images_names[uind].get_copymark().move_to(images_names);
                  dist.move_to(images,~0U);
                } else dist.move_to(images[uind].assign());
              }
          } else arg_error("distance");
          is_change = true;
          ++position;
          continue;
        }

        // Dilate.
        if (!std::strcmp("dilate",command)) {
          gmic_substitute_args(true);
          double sx = 3, sy = 3, sz = 1;
          unsigned int is_real = 0;
          boundary = 1;
          sep = 0;
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                            gmic_use_indices,&sep,&end)==2 && sep==']') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u%c",
                           indices,&boundary,&end)==2 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u%c",
                           indices,&boundary,&is_real,&end)==3) &&
              (ind=selection2cimg(indices,images.size(),images_names,"dilate")).height()==1 &&
              boundary<=3) {
            print(images,0,"Dilate image%s with kernel [%u] and %s boundary conditions, "
                  "in %s mode.",
                  gmic_selection.data(),
                  *ind,
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror",
                  is_real?"real":"binary");
            const CImg<T> kernel = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply(dilate(kernel,boundary,(bool)is_real),false);
          } else if ((cimg_sscanf(argument,"%lf%c",
                                  &sx,&end)==1) &&
                     sx>=0) {
            sx = cimg::round(sx);
            print(images,0,"Dilate image%s with kernel of size %g and neumann boundary conditions.",
                  gmic_selection.data(),
                  sx);
            cimg_forY(selection,l) gmic_apply(dilate((unsigned int)sx),true);
          } else if ((cimg_sscanf(argument,"%lf,%lf%c",
                                  &sx,&sy,&end)==2 ||
                      cimg_sscanf(argument,"%lf,%lf,%lf%c",
                                  &sx,&sy,&sz,&end)==3) &&
                     sx>=0 && sy>=0 && sz>=0) {
            sx = cimg::round(sx);
            sy = cimg::round(sy);
            sz = cimg::round(sz);
            print(images,0,"Dilate image%s with %gx%gx%g kernel and neumann boundary conditions.",
                  gmic_selection.data(),
                  sx,sy,sz);
            cimg_forY(selection,l) gmic_apply(dilate((unsigned int)sx,(unsigned int)sy,(unsigned int)sz),true);
          } else arg_error("dilate");
          is_change = true;
          ++position;
          continue;
        }

        // Patch-based smoothing.
        if (!std::strcmp("denoise",command)) {
          gmic_substitute_args(true);
          float sigma_s = 10, sigma_r = 10, smoothness = 1;
          unsigned int is_fast_approximation = 0;
          double psize = 5, rsize = 6;
          sep0 = sep1 = *argx = *argy = 0;
          if ((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-]%c",
                           gmic_use_indices,gmic_use_argx,&end)==2 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           indices,argx,gmic_use_argy,&end)==3 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],%lf%c",
                           indices,argx,argy,&psize,&end)==4 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],%lf,%lf%c",
                           indices,argx,argy,&psize,&rsize,&end)==5 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],%lf,%lf,%f%c",
                           indices,argx,argy,&psize,&rsize,&smoothness,&end)==6 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],%lf,%lf,%f,%u%c",
                           indices,argx,argy,&psize,&rsize,&smoothness,
                           &is_fast_approximation,&end)==7) &&
              (cimg_sscanf(argx,"%f%c",&sigma_s,&end)==1 ||
               (cimg_sscanf(argx,"%f%c%c",&sigma_s,&sep0,&end)==2 && sep0=='%')) &&
              (cimg_sscanf(argy,"%f%c",&sigma_r,&end)==1 ||
               (cimg_sscanf(argy,"%f%c%c",&sigma_r,&sep1,&end)==2 && sep1=='%')) &&
              (ind=selection2cimg(indices,images.size(),images_names,"denoise")).height()==1 &&
              sigma_s>=0 && sigma_r>=0 && psize>=0 && rsize>=0 && is_fast_approximation<=1) {
            psize = cimg::round(psize);
            rsize = cimg::round(rsize);
            print(images,0,"Denoise image%s using guide image [%u], with %gx%g patches, "
                  "standard deviations (%g%s,%g%s), lookup size %g and smoothness %g.",
                  gmic_selection.data(),*ind,psize,psize,
                  sigma_s,sep0=='%'?"%":"",
                  sigma_r,sep1=='%'?"%":"",
                  rsize,smoothness);
            const CImg<T> guide = gmic_image_arg(*ind);
            if (sep0=='%') sigma_s = -sigma_s;
            if (sep1=='%') sigma_r = -sigma_r;
            cimg_forY(selection,l)
              gmic_apply(blur_patch(guide,sigma_s,sigma_r,(unsigned int)psize,(unsigned int)rsize,smoothness,
                                    (bool)is_fast_approximation),false);
          } else if ((cimg_sscanf(argument,"%255[0-9.eE%+-]%c",
                                  gmic_use_argx,&end)==1 ||
                      cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                                  argx,gmic_use_argy,&end)==2 ||
                      cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%lf%c",
                                  argx,argy,&psize,&end)==3 ||
                      cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%lf,%lf%c",
                                  argx,argy,&psize,&rsize,&end)==4 ||
                      cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%lf,%lf,%f%c",
                                  argx,argy,&psize,&rsize,&smoothness,&end)==5 ||
                      cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%lf,%lf,%f,%u%c",
                                  argx,argy,&psize,&rsize,&smoothness,
                                  &is_fast_approximation,&end)==6) &&
                     (cimg_sscanf(argx,"%f%c",&sigma_s,&end)==1 ||
                      (cimg_sscanf(argx,"%f%c%c",&sigma_s,&sep0,&end)==2 && sep0=='%')) &&
                     (cimg_sscanf(argy,"%f%c",&sigma_r,&end)==1 ||
                      (cimg_sscanf(argy,"%f%c%c",&sigma_r,&sep1,&end)==2 && sep1=='%')) &&
                     sigma_s>=0 && sigma_r>=0 && psize>=0 && rsize>=0 && is_fast_approximation<=1) {
            psize = cimg::round(psize);
            rsize = cimg::round(rsize);
            print(images,0,"Denoise image%s with %gx%g patches, standard deviations (%g%s,%g%s), "
                  "lookup size %g and smoothness %g.",
                  gmic_selection.data(),psize,psize,
                  sigma_s,sep0=='%'?"%":"",
                  sigma_r,sep1=='%'?"%":"",
                  rsize,smoothness);
            if (sep0=='%') sigma_s = -sigma_s;
            if (sep1=='%') sigma_r = -sigma_r;
            cimg_forY(selection,l)
              gmic_apply(blur_patch(sigma_s,sigma_r,(unsigned int)psize,(unsigned int)rsize,smoothness,
                                    (bool)is_fast_approximation),false);
          } else arg_error("denoise");
          is_change = true;
          ++position;
          continue;
        }

        // Deriche filter.
        if (!std::strcmp("deriche",command)) {
          gmic_substitute_args(false);
          unsigned int order = 0;
          float sigma = 0;
          axis = sep = 0;
          boundary = 1;
          if ((cimg_sscanf(argument,"%f,%u,%c%c",&sigma,&order,&axis,&end)==3 ||
               (cimg_sscanf(argument,"%f%c,%u,%c%c",&sigma,&sep,&order,&axis,&end)==4 &&
                sep=='%') ||
               cimg_sscanf(argument,"%f,%u,%c,%u%c",&sigma,&order,&axis,&boundary,&end)==4 ||
               (cimg_sscanf(argument,"%f%c,%u,%c,%u%c",
                            &sigma,&sep,&order,&axis,&boundary,&end)==5 && sep=='%')) &&
              sigma>=0 && order<=2 && is_xyzc(axis) && boundary<=3) {
            print(images,0,"Apply %u-order Deriche filter on image%s, along axis '%c' with standard "
                  "deviation %g%s and %s boundary conditions.",
                  order,gmic_selection.data(),axis,
                  sigma,sep=='%'?"%":"",
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror");
            if (sep=='%') sigma = -sigma;
            cimg_forY(selection,l) gmic_apply(deriche(sigma,order,axis,boundary),true);
          } else arg_error("deriche");
          is_change = true;
          ++position;
          continue;
        }

        // Dijkstra algorithm.
        if (!std::strcmp("dijkstra",command)) {
          gmic_substitute_args(false);
          double snode = 0, enode = 0;
          if (cimg_sscanf(argument,"%lf,%lf%c",&snode,&enode,&end)==2 &&
              snode>=0 && enode>=0) {
            snode = cimg::round(snode);
            enode = cimg::round(enode);
            print(images,0,"Compute minimal path from adjacency matri%s%s with the "
                  "Dijkstra algorithm.",
                  selection.height()>1?"ce":"x",gmic_selection.data());
            unsigned int off = 0;
            cimg_forY(selection,l) {
              uind = selection[l] + off;
              CImg<T> path;
              if (is_get) {
                CImg<T> dist = gmic_check(images[uind]).get_dijkstra((unsigned int)snode,
                                                                     (unsigned int)enode,
                                                                     path);
                dist.move_to(images);
                path.move_to(images);
                images_names[uind].get_copymark().move_to(images_names);
                images_names.back().get_copymark().move_to(images_names);
              } else {
                gmic_check(images[uind]).dijkstra((unsigned int)snode,(unsigned int)enode,path);
                path.move_to(images,uind + 1);
                images_names[uind].get_copymark().move_to(images_names,uind + 1);
                ++off;
              }
            }
          } else arg_error("dijkstra");
          is_change = true;
          ++position;
          continue;
        }

        // Estimate displacement field.
        if (!std::strcmp("displacement",command)) {
          gmic_substitute_args(true);
          double nb_scales = 0, nb_iterations = 10000;
          float smoothness = 0.1f, precision = 5.f;
          unsigned int is_backward = 1;
          sep = *argx = 0;
          ind0.assign();
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                            gmic_use_indices,&sep,&end)==2 && sep==']') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%f%c",
                           indices,&smoothness,&end)==2 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%f,%f%c",
                           indices,&smoothness,&precision,&end)==3 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%f,%f,%lf%c",
                           indices,&smoothness,&precision,&nb_scales,&end)==4 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%f,%f,%lf,%lf%c",
                           indices,&smoothness,&precision,&nb_scales,&nb_iterations,&end)==5 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%f,%f,%lf,%lf,%u%c",
                           indices,&smoothness,&precision,&nb_scales,&nb_iterations,
                           &is_backward,&end)==6 ||
               (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%f,%f,%lf,%lf,%u,[%255[a-zA-Z0-9_.%+-]%c%c",
                            indices,&smoothness,&precision,&nb_scales,&nb_iterations,
                            &is_backward,gmic_use_argx,&sep,&end)==8 && sep==']')) &&
              (ind=selection2cimg(indices,images.size(),images_names,"displacement")).height()==1 &&
              precision>=0 && nb_scales>=0 && nb_iterations>=0 && is_backward<=1 &&
              (!*argx || (ind0=selection2cimg(argx,images.size(),images_names,"displacement")).height()==1)) {
            nb_scales = cimg::round(nb_scales);
            nb_iterations = cimg::round(nb_iterations);
            if (nb_scales) cimg_snprintf(argx,_argx.width(),"%g ",nb_scales); else std::strcpy(argx,"auto-");
            if (ind0) { gmic_use_argy; cimg_snprintf(argy,_argy.width()," with guide [%u]",*ind0); } else *argy = 0;

            print(images,0,"Estimate displacement field from source [%u] to image%s, with "
                  "%s smoothness %g, precision %g, %sscales, %g iteration%s, in %s direction%s.",
                  *ind,
                  gmic_selection.data(),
                  smoothness>=0?"isotropic":"anisotropic",cimg::abs(smoothness),
                  precision,
                  argx,
                  nb_iterations,nb_iterations!=1?"s":"",
                  is_backward?"backward":"forward",
                  argy);
            const CImg<T> source = gmic_image_arg(*ind);
            const CImg<T> constraints = ind0?gmic_image_arg(*ind0):CImg<T>::empty();
            cimg_forY(selection,l) gmic_apply(displacement(source,smoothness,precision,(unsigned int)nb_scales,
                                                           (unsigned int)nb_iterations,(bool)is_backward,
                                                           constraints),false);
          } else arg_error("displacement");
          is_change = true;
          ++position;
          continue;
        }

        // Delete file(s).
        if (!is_get && !std::strcmp("delete",item)) {
          gmic_substitute_args(false);
          const CImg<T> arg = CImg<char>::string(argument);
          const unsigned int pend = (unsigned int)arg.size();
          g_list_c.assign();
          for (unsigned int p = 0; p<pend; ) { // Retrieve list of filenames
            unsigned int np = p;
            while (np<pend && arg[np] && arg[np]!=',') ++np;
            if (np<pend) {
              CImg<T>(arg.data(p),1,++np - p,1,1,true).move_to(g_list_c);
              g_list_c.back().back() = 0;
            }
            p = np;
          }
          print(images,0,"Delete file%s '%s' (%u file%s).",
                g_list_c.size()!=1?"s":"",gmic_argument_text_printed(),
                g_list_c.size(),g_list_c.size()!=1?"s":"");
          cimglist_for(g_list_c,l) { strreplace_fw(g_list_c[l]); std::remove(g_list_c[l]); }
          g_list_c.assign();
          ++position;
          continue;
        }

        // Display.
        if (!is_get && !std::strcmp("display",command)) {
          gmic_substitute_args(false);
          *argx = *argy = *argz = sep = sep0 = sep1 = 0;
          value = value0 = value1 = 0;
          exit_on_anykey = 0;
          if (((cimg_sscanf(argument,"%255[0-9.eE%+-]%c",
                            gmic_use_argx,&end)==1) ||
               (cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                            argx,gmic_use_argy,&end)==2) ||
               (cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                            argx,argy,gmic_use_argz,&end)==3) ||
               (cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%u%c",
                            argx,argy,argz,&exit_on_anykey,&end)==4)) &&
              (cimg_sscanf(argx,"%lf%c",&value,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&value,&sep,&end)==2 && sep=='%')) &&
              (!*argy ||
               cimg_sscanf(argy,"%lf%c",&value0,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&value0,&sep0,&end)==2 && sep0=='%')) &&
              (!*argz ||
               cimg_sscanf(argz,"%lf%c",&value1,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&value1,&sep1,&end)==2 && sep1=='%')) &&
              value>=0 && value0>=0 && value1>=0 && exit_on_anykey<=1) {
            if (!*argx) { value = 50; sep = '%'; }
            if (!*argy) { value0 = 50; sep0 = '%'; }
            if (!*argz) { value1 = 50; sep1 = '%'; }
            ++position;
          } else { value = value0 = value1 = 50; sep = sep0 = sep1 = '%'; }

          unsigned int XYZ[3];
          if (selection.height()>=1) {
            CImg<T> &img = images[selection[0]];
            XYZ[0] = (unsigned int)cimg::cut(cimg::round(sep=='%'?(img.width() - 1)*value/100:value),
                                             0.,img.width() - 1.);
            XYZ[1] = (unsigned int)cimg::cut(cimg::round(sep0=='%'?(img.height() - 1)*value0/100:value0),
                                             0.,img.height() - 1.);
            XYZ[2] = (unsigned int)cimg::cut(cimg::round(sep1=='%'?(img.depth() - 1)*value1/100:value1),
                                             0.,img.depth() - 1.);
          }
          ++verbosity;
          display_images(images,images_names,selection,XYZ,exit_on_anykey);
          --verbosity;
          is_change = false;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'e...'
        //-----------------------------
      gmic_commands_e :
        if (check_elif && !std::strcmp("elif",item)) // Redirect 'elif'
          goto gmic_commands_others;

        // Else and eluded elif.
        if (!is_get && (!std::strcmp("else",item) || (!check_elif && !std::strcmp("elif",item)))) {
          const CImg<char> &s = callstack.back();
          if (s[0]!='*' || s[1]!='i')
            error(true,images,0,0,
                  "Command '%s': Not associated to a 'if' command within the same scope.",
                  item);
          check_elif = false;
          if (is_very_verbose) print(images,0,"Reach 'else' block.");
          for (int nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
            it = commands_line[position];
            if (*it==1)
              is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
            else {
              it+=*it=='-';
              if (!std::strcmp("if",it)) ++nb_levels;
              else if (!std::strcmp("fi",it)) { if (!--nb_levels) --position; }
            }
          }
          continue;
        }

        // Evaluate expression.
        if (!std::strcmp("eval",command)) {
          if (is_get && !is_selection)
            error(true,images,0,0,
                  "Command 'eval': Image selection is missing.");
          gmic_substitute_args(false);
          gmic_argument_text_printed();
          if (*_argument_text=='\'') cimg::strpare(_argument_text,'\'',true,false);
          name.assign(argument,(unsigned int)std::strlen(argument) + 1);
          cimg::strpare(name,'\'',true,false);
          strreplace_fw(name);

          if (!is_selection) { // No selection -> single evaluation
            print(images,0,"Evaluate expression '%s' and assign it to status.",
                  _argument_text.data());
            CImg<T> &img = images.size()?images.back():CImg<T>::empty();
            CImg<double> output;
            img.eval(output,name,0,0,0,0,&images);
            if (output.height()>1) // Vector-valued result
              output.value_string().move_to(status);
            else { // Scalar result
              gmic_use_formula;
              cimg_snprintf(formula,_formula.width(),"%.17g",*output);
              CImg<char>::string(formula).move_to(status);
            }
          } else { // Selection -> loop over images
            print(images,0,"Evaluate expression '%s' looped over image%s.",
                  _argument_text.data(),
                  gmic_selection.data());
            cimg_forY(selection,l) gmic_apply(gmic_eval(name.data(),images),false);
            is_change = true;
          }
          ++position;
          continue;
        }

        // Echo.
        if (is_command_echo) {
          if (verbosity>=0 || is_debug || is_get) {
            gmic_substitute_args(false);
            name.assign(argument,(unsigned int)std::strlen(argument) + 1);
            cimg::strunescape(name);
            const int _verbosity = ++verbosity;
            std::FILE *_file = 0;
            if (is_get) { _file = cimg::output(); verbosity = 1; cimg::output(stdout); }
            if (is_selection) print(images,&selection,"%s",name.data());
            else print(images,&CImg<unsigned int>::empty(),"%s",name.data());
            if (is_get) { verbosity = _verbosity; cimg::output(_file); }
            --verbosity;
          }
          ++position;
          continue;
        }

        // Exec.
        if (!is_get && !std::strcmp("exec",item)) {
          gmic_substitute_args(false);
          name.assign(argument,(unsigned int)std::strlen(argument) + 1);
          const char *arg_exec_text = gmic_argument_text_printed();
          char *arg_exec = name;
          cimg::strunescape(arg_exec);
          strreplace_fw(arg_exec);

          is_verbose = true; // is_verbose
          if ((*arg_exec=='0' || *arg_exec=='1') && arg_exec[1]==',') {
            is_verbose = (*arg_exec=='1');
            arg_exec+=2; arg_exec_text+=2;
          }

#ifdef gmic_noexec
          print(images,0,"Execute external command '%s' %s (skipped, no exec allowed).",
                arg_exec_text,
                is_verbose?"in verbose mode":"");
#else // #ifdef gmic_noexec
          print(images,0,"Execute external command '%s' %s",
                arg_exec_text,
                is_verbose?"in verbose mode.\n":".");
          cimg::mutex(31);
          const int errcode = cimg::system(arg_exec,0,is_verbose);
          cimg::mutex(31,0);
          gmic_use_title;
          cimg_snprintf(title,_title.width(),"%d",errcode);
          CImg<char>::string(title).move_to(status);
          if (errcode) print(images,0,"Command 'exec' returned error code %d.",
                             errcode);
#endif // #ifdef gmic_noexec
          ++position;
          continue;
        }

        // Error.
        if (is_command_error) {
          gmic_substitute_args(false);
          name.assign(argument,(unsigned int)std::strlen(argument) + 1);
          cimg::strunescape(name);
          if (is_selection) error(true,images,&selection,0,"%s",name.data());
          else error(true,images,&CImg<unsigned int>::empty(),0,"%s",name.data());
        }

        // Invert endianness.
        if (!std::strcmp("endian",command)) {
          gmic_substitute_args(false);
          if (!std::strcmp(argument,"bool") ||
              !std::strcmp(argument,"uint8") || !std::strcmp(argument,"int8") ||
              !std::strcmp(argument,"uint16") || !std::strcmp(argument,"int16") ||
              !std::strcmp(argument,"uint32") || !std::strcmp(argument,"int32") ||
              !std::strcmp(argument,"uint64") || !std::strcmp(argument,"int64") ||
              !std::strcmp(argument,"float32") || !std::strcmp(argument,"float64")) {
            print(images,0,"Invert data endianness of image%s, with assumed pixel type '%s'.",
                  gmic_selection.data(),argument);
            ++position;
          } else print(images,0,"Invert data endianness of image%s.",
                       gmic_selection.data());
          cimg_forY(selection,l) gmic_apply(gmic_invert_endianness(argument),true);
          is_change = true;
          continue;
        }

        // Exponential.
        gmic_simple_command("exp",exp,"Compute pointwise exponential of image%s.");

        // Error function.
        gmic_simple_command("erf",erf,"Compute pointwise error function of image%s.");

        // Test equality.
        gmic_arithmetic_command("eq",
                                operator_eq,
                                "Compute boolean equality between image%s and %g%s",
                                gmic_selection.data(),value,ssep,T,
                                operator_eq,
                                "Compute boolean equality between image%s and image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_eq,
                                "Compute boolean equality between image%s and expression %s'",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute boolean equality between image%s");

        // Draw ellipse.
        if (!std::strcmp("ellipse",command)) {
          gmic_substitute_args(false);
          double x = 0, y = 0, R = 0, r = 0, angle = 0;
          sep = sepx = sepy = sepz = sepc = *argx = *argy = *argz = *argc = *color = 0;
          pattern = ~0U; opacity = 1;
          if ((cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           gmic_use_argx,gmic_use_argy,gmic_use_argz,&end)==3 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-]%c",
                           argx,argy,argz,gmic_use_argc,&end)==4 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-],%lf%c",
                           argx,argy,argz,argc,&angle,&end)==5 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-],%lf,%f%c",
                           argx,argy,argz,argc,&angle,&opacity,&end)==6 ||
               (cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                            "%255[0-9.eE%+-],%lf,%f,0%c%x%c",
                            argx,argy,argz,argc,&angle,&opacity,&sep,&pattern,
                            &end)==8 &&
                sep=='x') ||
               (cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                            "%255[0-9.eE%+-],%lf,%f,0%c%x,%4095[0-9.eEinfa,+-]%c",
                            argx,argy,argz,argc,&angle,&opacity,&sep,
                            &pattern,gmic_use_color,&end)==9 &&
                sep=='x') ||
               (cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                            "%255[0-9.eE%+-],%lf,%f,%4095[0-9.eEinfa,+-]%c",
                            argx,argy,argz,argc,&angle,&opacity,color,&end)==7)) &&
              (cimg_sscanf(argx,"%lf%c",&x,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (cimg_sscanf(argy,"%lf%c",&y,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
              (cimg_sscanf(argz,"%lf%c",&R,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&R,&sepz,&end)==2 && sepz=='%')) &&
              (!*argc ||
               cimg_sscanf(argc,"%lf%c",&r,&end)==1 ||
               (cimg_sscanf(argc,"%lf%c%c",&r,&sepc,&end)==2 && sepc=='%'))) {
            if (!*argc) r = R;
            print(images,0,"Draw %s ellipse at (%g%s,%g%s) with radii (%g%s,%g%s) on image%s, "
                  "with orientation %g deg., opacity %g and color (%s).",
                  sep=='x'?"outlined":"filled",
                  x,sepx=='%'?"%":"",
                  y,sepy=='%'?"%":"",
                  R,sepz=='%'?"%":"",
                  r,sepc=='%'?"%":"",
                  gmic_selection.data(),
                  angle,
                  opacity,
                  *color?color:"default");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              g_img.assign(img.spectrum(),1,1,1,(T)0).fill_from_values(color,true);
              const float rmax = std::sqrt((float)cimg::sqr(img.width()) +
                                           cimg::sqr(img.height()));
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.width() - 1)/100:x),
                ny = (int)cimg::round(sepy=='%'?y*(img.height() - 1)/100:y);
              const float
                nR = (float)(sepz=='%'?R*rmax/100:R),
                nr = (float)(sepc=='%'?r*rmax/100:r);
              if (sep=='x') {
                gmic_apply(draw_ellipse(nx,ny,nR,nr,(float)angle,g_img.data(),opacity,pattern),true);
              } else {
                gmic_apply(draw_ellipse(nx,ny,nR,nr,(float)angle,g_img.data(),opacity),true);
              }
            }
          } else arg_error("ellipse");
          g_img.assign();
          is_change = true;
          ++position;
          continue;
        }

        // Equalize.
        if (!std::strcmp("equalize",command)) {
          gmic_substitute_args(false);
          double nb_levels = 256;
          bool no_min_max = false;
          sep = sep0 = sep1 = 0;
          value0 = value1 = 0;
          if (((cimg_sscanf(argument,"%lf%c",
                            &nb_levels,&end)==1 && (no_min_max=true)) ||
               ((cimg_sscanf(argument,"%lf%c%c",
                             &nb_levels,&sep,&end)==2 && sep=='%') && (no_min_max=true)) ||
               cimg_sscanf(argument,"%lf,%lf,%lf%c",
                           &nb_levels,&value0,&value1,&end)==3 ||
               (cimg_sscanf(argument,"%lf%c,%lf,%lf%c",
                            &nb_levels,&sep,&value0,&value1,&end)==4 && sep=='%') ||
               (cimg_sscanf(argument,"%lf,%lf%c,%lf%c",
                            &nb_levels,&value0,&sep0,&value1,&end)==4 && sep0=='%') ||
               (cimg_sscanf(argument,"%lf%c,%lf%c,%lf%c",
                            &nb_levels,&sep,&value0,&sep0,&value1,&end)==5 && sep=='%' &&
                sep0=='%') ||
               (cimg_sscanf(argument,"%lf,%lf,%lf%c%c",
                            &nb_levels,&value0,&value1,&sep1,&end)==4 && sep1=='%') ||
               (cimg_sscanf(argument,"%lf%c,%lf,%lf%c%c",
                            &nb_levels,&sep,&value0,&value1,&sep1,&end)==5 && sep=='%' &&
                sep1=='%') ||
               (cimg_sscanf(argument,"%lf,%lf%c,%lf%c%c",
                            &nb_levels,&value0,&sep0,&value1,&sep1,&end)==5 && sep0=='%' &&
                sep1=='%') ||
               (cimg_sscanf(argument,"%lf%c,%lf%c,%lf%c%c",
                            &nb_levels,&sep,&value0,&sep0,&value1,&sep1,&end)==6 && sep=='%' &&
                sep0=='%' && sep1=='%')) &&
              nb_levels>=0.5) { nb_levels = cimg::round(nb_levels); ++position; }
          else { nb_levels = 256; value0 = 0; value1 = 100; sep = 0; sep0 = sep1 = '%'; }
          if (no_min_max) { value0 = 0; value1 = 100; sep0 = sep1 = '%'; }
          print(images,0,"Equalize histogram of image%s, with %g%s levels in range [%g%s,%g%s].",
                gmic_selection.data(),
                nb_levels,sep=='%'?"%":"",
                value0,sep0=='%'?"%":"",
                value1,sep1=='%'?"%":"");
          cimg_forY(selection,l) {
            CImg<T>& img = gmic_check(images[selection[l]]);
            nvalue0 = value0; nvalue1 = value1;
            vmin = vmax = 0;
            if (sep0=='%' || sep1=='%') {
              if (img) vmax = (double)img.max_min(vmin);
              if (sep0=='%') nvalue0 = vmin + (vmax - vmin)*value0/100;
              if (sep1=='%') nvalue1 = vmin + (vmax - vmin)*value1/100;
            }
            const unsigned int
              _nb_levels = std::max(1U,
                                    (unsigned int)cimg::round(sep=='%'?
                                                              nb_levels*(1 + nvalue1 - nvalue0)/100:
                                                              nb_levels));
            gmic_apply(equalize(_nb_levels,(T)nvalue0,(T)nvalue1),true);
          }
          is_change = true;
          continue;
        }

        // Erode.
        if (!std::strcmp("erode",command)) {
          gmic_substitute_args(true);
          unsigned int is_real = 0;
          double sx = 3, sy = 3, sz = 1;
          boundary = 1;
          sep = 0;
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                            gmic_use_indices,&sep,&end)==2 && sep==']') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u%c",
                           indices,&boundary,&end)==2 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u%c",
                           indices,&boundary,&is_real,&end)==3) &&
              (ind=selection2cimg(indices,images.size(),images_names,"erode")).height()==1 &&
              boundary<=3) {
            print(images,0,"Erode image%s with kernel [%u] and %s boundary conditions, "
                  "in %s mode.",
                  gmic_selection.data(),
                  *ind,
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror",
                  is_real?"real":"binary");
            const CImg<T> kernel = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply(erode(kernel,boundary,(bool)is_real),false);
          } else if ((cimg_sscanf(argument,"%lf%c",
                                  &sx,&end)==1) &&
                     sx>=0) {
            sx = cimg::round(sx);
            print(images,0,"Erode image%s with kernel of size %g and neumann boundary conditions.",
                  gmic_selection.data(),
                  sx);
            cimg_forY(selection,l) gmic_apply(erode((unsigned int)sx),true);
          } else if ((cimg_sscanf(argument,"%lf,%lf%c",
                                  &sx,&sy,&end)==2 ||
                      cimg_sscanf(argument,"%lf,%lf,%lf%c",
                                  &sx,&sy,&sz,&end)==3) &&
                     sx>=0 && sy>=0 && sz>=0) {
            sx = cimg::round(sx);
            sy = cimg::round(sy);
            sz = cimg::round(sz);
            print(images,0,"Erode image%s with %gx%gx%g kernel and neumann boundary conditions.",
                  gmic_selection.data(),
                  sx,sy,sz);
            cimg_forY(selection,l) gmic_apply(erode((unsigned int)sx,(unsigned int)sy,(unsigned int)sz),true);
          } else arg_error("erode");
          is_change = true;
          ++position;
          continue;
        }

        // Eigenvalues/eigenvectors.
        if (!std::strcmp("eigen",command)) {
          print(images,0,"Compute eigen-values/vectors of symmetric matri%s or matrix field%s.",
                selection.height()>1?"ce":"x",gmic_selection.data());
          unsigned int off = 0;
          cimg_forY(selection,l) {
            uind = selection[l] + off;
            CImg<float> val, vec;
            gmic_check(images[uind]).gmic_symmetric_eigen(val,vec);
            if (is_get) {
              val.move_to(images);
              vec.move_to(images);
              images_names[uind].get_copymark().move_to(images_names);
              images_names.back().get_copymark().move_to(images_names);
            } else {
              val.move_to(images[uind].assign());
              vec.move_to(images,uind + 1);
              images_names[uind].get_copymark().move_to(images_names,uind + 1);
              ++off;
            }
          }
          is_change = true;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'f...'
        //-----------------------------
      gmic_commands_f :
        if (command[1]=='f' && command[2]=='t' && !command[3]) goto gmic_commands_e; // Redirect 'fft'

        // Fi.
        if (!is_get && !std::strcmp("fi",item)) {
          const CImg<char> &s = callstack.back();
          if (s[0]!='*' || s[1]!='i')
            error(true,images,0,0,
                  "Command 'endif': Not associated to a 'if' command within the same scope.");
          if (is_very_verbose) print(images,0,"End 'if...endif' block.");
          check_elif = false;
          callstack.remove();
          continue;
        }

        // For.
        if (!is_get && !std::strcmp("for",item)) {
          gmic_substitute_args(false);
          is_cond = check_cond(argument,images,"for");
          const bool is_first = !nb_fordones || fordones(0U,nb_fordones - 1)!=position_item;
          if (is_very_verbose)
            print(images,0,"%s %s -> condition '%s' %s.",
                  !is_first?"Go back to":is_cond?"Start":"Skip",
                  is_first?"'for...done' block":"'for' command",
                  gmic_argument_text_printed(),
                  is_cond?"holds":"does not hold");
          ++position;

          if (!is_cond) {
            int nb_levels = 0;
            for (nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
              it = commands_line[position];
              if (*it==1)
                is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
              else {
                _is_get = *it=='+';
                it+=(_is_get || *it=='-');
                gmic_if_flr ++nb_levels; gmic_elif_flr --nb_levels;
              }
            }
            if (nb_levels)
              error(true,images,0,0,
                    "Command 'for': Missing associated 'done' command.");
            if (!is_first) { --nb_fordones; callstack.remove(); }
          } else if (is_first) {
            if (is_debug_info && debug_line!=~0U) {
              gmic_use_argx;
              cimg_snprintf(argx,_argx.width(),"*for#%u",debug_line);
              CImg<char>::string(argx).move_to(callstack);
            } else CImg<char>::string("*for").move_to(callstack);
            if (nb_fordones>=fordones._height) fordones.resize(3,std::max(2*fordones._height,8U),1,1,0);
            unsigned int *const fd = fordones.data(0,nb_fordones++);
            fd[0] = position_item;
            fd[1] = 0;
            fd[2] = debug_line;
          }
          is_lbrace_command = true;
          continue;
        }

        // Foreach.
        if (!std::strcmp("foreach",command)) {
          if (!selection) {
            if (is_very_verbose) print(images,0,"Skip 'foreach...done' block.");
            int nb_levels = 0;
            for (nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
              it = commands_line[position];
              if (*it==1)
                is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
              else {
                _is_get = *it=='+';
                it+=(_is_get || *it=='-');
                gmic_if_flr ++nb_levels; gmic_elif_flr --nb_levels;
              }
            }
            if (nb_levels)
              error(true,images,0,0,
                    "Command 'foreach': Missing associated 'done' command.");
          } else {
            if (is_debug_info && debug_line!=~0U) {
              gmic_use_argx;
              cimg_snprintf(argx,_argx.width(),"*foreach#%u",debug_line);
              CImg<char>::string(argx).move_to(callstack);
            } else CImg<char>::string("*foreach").move_to(callstack);
            if (is_very_verbose)
              print(images,0,"Start 'foreach...done' block, with image%s.",
                    gmic_selection.data());

            if (nb_foreachdones>=foreachdones._height)
              foreachdones.resize(3,std::max(2*foreachdones._height,8U),1,1,0);
            unsigned int *fed = foreachdones.data(0,nb_foreachdones++);
            fed[0] = 0; // Iteration counter
            fed[1] = selection._height;
            fed[2] = debug_line;

            const unsigned int _position = position;
            int off = 0;

            cimg_forY(selection,l) {
              uind = selection[l] + off;
              if (is_get) {
                g_list.assign(images[uind]);
                g_list_c.assign(images_names[uind]);
              } else {
                cimg::mutex(27);
                if (images[uind].is_shared()) g_list.assign(images[uind]);
                else
                  if ((images[uind].width() || images[uind].height()) && !images[uind]._spectrum) {
                    selection2string(selection,images_names,1,name);
                    error(true,images,0,0,
                          "Command 'foreach': Invalid selection%s "
                          "(image [%u] is already used in another thread).",
                          name.data() + (*name=='s'?1:0),uind);
                  }
                images[uind].move_to(g_list);
                g_list_c.assign(images_names[uind]);

                // Small hack to be able to track images of the selection passed to the new environment.
                std::memcpy(&images[uind]._width,&g_list[0]._data,sizeof(void*));
                images[uind]._spectrum = 0;
                cimg::mutex(27,0);
              }

              is_lbrace_command = true;
              gmic_exception exception;
              try {
                if (next_debug_line!=~0U) { debug_line = next_debug_line; next_debug_line = ~0U; }
                if (next_debug_filename!=~0U) { debug_filename = next_debug_filename; next_debug_filename = ~0U; }
                _run(commands_line,position = _position,g_list,g_list_c,images,images_names,variables_sizes,is_noarg,0,
                     command_selection,false);
              } catch (gmic_exception &e) {
                check_elif = false;
                int nb_levels = 0;
                for (nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
                  it = commands_line[position];
                  if (*it==1)
                    is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
                  else {
                    _is_get = *it=='+';
                    if (*it==1)
                      is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
                    else {
                      it+=(_is_get || *it=='-');
                      gmic_if_flr ++nb_levels; gmic_elif_flr --nb_levels;
                    }
                  }
                }
                cimg::swap(exception._command,e._command);
                cimg::swap(exception._message,e._message);
              }

              // Transfer back images to saved list of images.
              cimg::mutex(27);
              if (is_get) {
                g_list.move_to(images,~0U);
                cimglist_for(g_list_c,i) g_list_c[i].get_copymark().move_to(g_list_c[i]);
                g_list_c.move_to(images_names,~0U);
              } else {
                if (!g_list) {
                  images.remove(uind);
                  images_names.remove(uind);
                  --off;
                } else if (g_list.size()==1) {
                  g_list[0].move_to(images[uind]);
                  g_list_c[0].move_to(images_names[uind]);
                } else {
                  images.remove(uind);
                  images_names.remove(uind);
                  off+=g_list.width() - 1;
                  g_list.move_to(images,uind);
                  g_list_c.move_to(images_names,uind);
                }
              }
              g_list.assign();
              g_list_c.assign();

              if (!exception._message && nb_foreachdones>0 && foreachdones._height>=nb_foreachdones) {
                fed = foreachdones.data(0,nb_foreachdones - 1);
                if (!fed[1]) l = selection.height(); // when case break() happened, force loop to stop
                else { ++fed[0]; --fed[1]; }
                next_debug_line = fed[2];
                next_debug_filename = debug_filename;
              }

              cimg::mutex(27,0);
              if (exception._message) throw exception;
            }

            --nb_foreachdones;
            callstack.remove();
          }
          continue;
        }

        // Fill.
        if (!std::strcmp("fill",command)) {
          gmic_substitute_args(true);
          sep = *indices = 0;
          value = 0;
          if (cimg_sscanf(argument,"%lf%c",
                          &value,&end)==1) {
            print(images,0,"Fill image%s with %g.",
                  gmic_selection.data(),
                  value);
            cimg_forY(selection,l) gmic_apply(fill((T)value),false);
          } else if (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep,&end)==2 &&
                     sep==']' &&
                     (ind=selection2cimg(indices,images.size(),images_names,"fill")).height()==1) {
            print(images,0,"Fill image%s with values from image [%u].",
                  gmic_selection.data(),
                  *ind);
            const CImg<T> values = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply(fill(values),false);
          } else {
            gmic_argument_text_printed();
            if (*_argument_text=='\'') cimg::strpare(_argument_text,'\'',true,false);
            print(images,0,"Fill image%s with expression '%s'.",
                  gmic_selection.data(),
                  _argument_text.data());
            name.assign(argument,(unsigned int)std::strlen(argument) + 1);
            cimg::strpare(name,'\'',true,false);
            strreplace_fw(name);
            cimg_forY(selection,l) gmic_apply(gmic_fill(name.data(),images),false);
          }
          is_change = true;
          ++position;
          continue;
        }

        // Flood fill.
        if (!std::strcmp("flood",command)) {
          gmic_substitute_args(false);
          double x = 0, y = 0, z = 0;
          float tolerance = 0;
          sepx = sepy = sepz = *argx = *argy = *argz = *color = 0;
          is_high_connectivity = 0;
          opacity = 1;
          if ((cimg_sscanf(argument,"%255[0-9.eE%+-]%c",
                           gmic_use_argx,&end)==1 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           argx,gmic_use_argy,&end)==2 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           argx,argy,gmic_use_argz,&end)==3 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eEinfa%+-],%f%c",
                           argx,argy,argz,&tolerance,&end)==4 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eEinfa%+-],%f,%u%c",
                           argx,argy,argz,&tolerance,&is_high_connectivity,&end)==5 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eEinfa%+-],%f,%u,%f%c",
                           argx,argy,argz,&tolerance,&is_high_connectivity,&opacity,&end)==6 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eEinfa%+-],%f,%u,%f,"
                           "%4095[0-9.eEinfa,+-]%c",
                           argx,argy,argz,&tolerance,&is_high_connectivity,
                           &opacity,gmic_use_color,&end)==7) &&
              (cimg_sscanf(argx,"%lf%c",&x,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy ||
               cimg_sscanf(argy,"%lf%c",&y,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
              (!*argz ||
               cimg_sscanf(argz,"%lf%c",&z,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&z,&sepz,&end)==2 && sepz=='%')) &&
              tolerance>=0) {
            print(images,0,
                  "Flood fill image%s from (%g%s,%g%s,%g%s), with tolerance %g, %s connectivity, "
                  "opacity %g and color (%s).",
                  gmic_selection.data(),
                  x,sepx=='%'?"%":"",
                  y,sepy=='%'?"%":"",
                  z,sepz=='%'?"%":"",
                  tolerance,
                  is_high_connectivity?"high":"low",
                  opacity,
                  *color?color:"default");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              g_img.assign(img.spectrum(),1,1,1,(T)0).fill_from_values(color,true);
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.width() - 1)/100:x),
                ny = (int)cimg::round(sepy=='%'?y*(img.height() - 1)/100:y),
                nz = (int)cimg::round(sepz=='%'?z*(img.depth() - 1)/100:z);
              gmic_apply(draw_fill(nx,ny,nz,g_img.data(),opacity,tolerance,(bool)is_high_connectivity),false);
            }
          } else arg_error("flood");
          g_img.assign();
          is_change = true;
          ++position;
          continue;
        }

        // List of directory files.
        if (!is_get && !std::strcmp("files",item)) {
          gmic_substitute_args(false);
          unsigned int mode = 5;
          if ((*argument>='0' && *argument<='5') &&
              argument[1]==',' && argument[2]) {
            mode = (unsigned int)(*argument - '0');
            argument+=2;
          }
          const unsigned int _mode = mode%3;
          print(images,0,"Get list of %s from location '%s'.",
                _mode==0?"files":_mode==1?"folders":"files and folders",
                argument);
          g_list_c = cimg::files(argument,true,_mode,mode>=3);
          cimglist_for(g_list_c,l) {
            strreplace_bw(g_list_c[l]);
            g_list_c[l].back() = ',';
          }
          if (g_list_c) {
            g_list_c.back().back() = 0;
            (g_list_c>'x').move_to(status);
          } else status.assign();
          g_list_c.assign();
          ++position;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'g...'
        //-----------------------------
      gmic_commands_g :

        // Greater or equal.
        gmic_arithmetic_command("ge",
                                operator_ge,
                                "Compute boolean 'greater or equal than' between image%s and %g%s",
                                gmic_selection.data(),value,ssep,T,
                                operator_ge,
                                "Compute boolean 'greater or equal than' between image%s "
                                "and image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_ge,
                                "Compute boolean 'greater or equal than' between image%s "
                                "and expression %s'",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute boolean 'greater or equal than' between image%s");

        // Greater than.
        gmic_arithmetic_command("gt",
                                operator_gt,
                                "Compute boolean 'greater than' between image%s and %g%s",
                                gmic_selection.data(),value,ssep,T,
                                operator_gt,
                                "Compute boolean 'greater than' between image%s and image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_gt,
                                "Compute boolean 'greater than' between image%s and expression %s'",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute boolean 'greater than' between image%s");

        // Guided filter.
        if (!std::strcmp("guided",command)) {
          gmic_substitute_args(true);
          float radius = 0, regularization = 0;
          sep0 = sep1 = 0;
          if (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                          gmic_use_indices,gmic_use_argx,gmic_use_argy,&end)==3 &&
              (cimg_sscanf(argx,"%f%c",&radius,&end)==1 ||
               (cimg_sscanf(argx,"%f%c%c",&radius,&sep0,&end)==2 && sep0=='%')) &&
              (cimg_sscanf(argy,"%f%c",&regularization,&end)==1 ||
               (cimg_sscanf(argy,"%f%c%c",&regularization,&sep1,&end)==2 && sep1=='%')) &&
              (ind=selection2cimg(indices,images.size(),images_names,"guided")).height()==1 &&
              radius>=0 && regularization>=0) {
            print(images,0,"Apply guided filter on image%s, with guide image [%u], "
                  "radius %g%s and regularization %g%s.",
                  gmic_selection.data(),
                  *ind,
                  radius,sep0=='%'?"%":"",
                  regularization,sep1=='%'?"%":"");
            const CImg<T> guide = gmic_image_arg(*ind);
            if (sep0=='%') radius = -radius;
            if (sep1=='%') regularization = -regularization;
            cimg_forY(selection,l) gmic_apply(blur_guided(guide,radius,regularization),false);
          } else if (cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                                 gmic_use_argx,gmic_use_argy,&end)==2 &&
                     (cimg_sscanf(argx,"%f%c",&radius,&end)==1 ||
                      (cimg_sscanf(argx,"%f%c%c",&radius,&sep0,&end)==2 && sep0=='%')) &&
                     (cimg_sscanf(argy,"%f%c",&regularization,&end)==1 ||
                      (cimg_sscanf(argy,"%f%c%c",&regularization,&sep1,&end)==2 && sep1=='%')) &&
                     radius>=0 && regularization>=0) {
            print(images,0,"Apply guided filter on image%s, with radius %g%s and regularization %g%s.",
                  gmic_selection.data(),
                  radius,sep0=='%'?"%":"",
                  regularization,sep1=='%'?"%":"");
            if (sep0=='%') radius = -radius;
            if (sep1=='%') regularization = -regularization;
            cimg_forY(selection,l) gmic_apply(blur_guided(images[selection[l]],radius,regularization),false);
          } else arg_error("guided");
          is_change = true;
          ++position;
          continue;
        }

        // Draw graph.
        if (!std::strcmp("graph",command)) {
          gmic_substitute_args(true);
          double ymin = 0, ymax = 0, xmin = 0, xmax = 0, resolution = 65536;
          unsigned int plot_type = 1, vertex_type = 1;
          *formula = *color = sep = sep1 = 0;
          pattern = ~0U; opacity = 1;
          if (((cimg_sscanf(argument,"'%1023[^']%c%c",
                            gmic_use_formula,&sep,&end)==2 && sep=='\'') ||
               cimg_sscanf(argument,"'%1023[^']',%lf%c",
                           formula,&resolution,&end)==2 ||
               cimg_sscanf(argument,"'%1023[^']',%lf,%u%c",
                           formula,&resolution,&plot_type,&end)==3 ||
               cimg_sscanf(argument,"'%1023[^']',%lf,%u,%u%c",
                           formula,&resolution,&plot_type,&vertex_type,&end)==4 ||
               cimg_sscanf(argument,"'%1023[^']',%lf,%u,%u,%lf,%lf%c",
                           formula,&resolution,&plot_type,&vertex_type,&xmin,&xmax,&end)==6 ||
               cimg_sscanf(argument,"'%1023[^']',%lf,%u,%u,%lf,%lf,%lf,%lf%c",
                           formula,&resolution,&plot_type,&vertex_type,&xmin,&xmax,
                           &ymin,&ymax,&end)==8 ||
               cimg_sscanf(argument,"'%1023[^']',%lf,%u,%u,%lf,%lf,%lf,%lf,%f%c",
                           formula,&resolution,&plot_type,&vertex_type,
                           &xmin,&xmax,&ymin,&ymax,&opacity,&end)==9 ||
               (cimg_sscanf(argument,"'%1023[^']',%lf,%u,%u,%lf,%lf,%lf,%lf,%f,0%c%x%c",
                            formula,&resolution,&plot_type,&vertex_type,&xmin,&xmax,
                            &ymin,&ymax,&opacity,&sep1,&pattern,&end)==11 && sep1=='x') ||
               (cimg_sscanf(argument,"'%1023[^']',%lf,%u,%u,%lf,%lf,%lf,%lf,%f,%4095[0-9.eEinfa,+-]%c",
                            formula,&resolution,&plot_type,&vertex_type,&xmin,&xmax,&ymin,&ymax,
                            &opacity,gmic_use_color,&end)==10 && (bool)(pattern=~0U)) ||
               (cimg_sscanf(argument,"'%1023[^']',%lf,%u,%u,%lf,%lf,%lf,%lf,%f,0%c%x,"
                            "%4095[0-9.eEinfa,+-]%c",
                            formula,&resolution,&plot_type,&vertex_type,&xmin,&xmax,
                            &ymin,&ymax,&opacity,&sep1,&pattern,&(*color=0),&end)==12 &&
                sep1=='x')) &&
              resolution>0 && plot_type<=3 && vertex_type<=7) {
            resolution = cimg::round(resolution);
            strreplace_fw(formula);
            print(images,0,
                  "Draw graph of formula '%s' on image%s, with resolution %g, %s contours, "
                  "%s vertices, x-range = (%g,%g), y-range = (%g,%g), opacity %g, "
                  "pattern 0x%x and color (%s).",
                  formula,
                  gmic_selection.data(),
                  resolution,
                  plot_type==0?"no":plot_type==1?"linear":plot_type==2?"spline":"bar",
                  vertex_type==0?"no":vertex_type==1?"dot":vertex_type==2?"straight cross":
                  vertex_type==3?"diagonal cross":vertex_type==4?"filled circle":
                  vertex_type==5?"outlined circle":vertex_type==6?"square":"diamond",
                  xmin,xmax,
                  ymin,ymax,
                  opacity,pattern,
                  *color?color:"default");
            if (xmin==0 && xmax==0) { xmin = -4; xmax = 4; }
            if (!plot_type && !vertex_type) plot_type = 1;
            if (resolution<1) resolution = 65536;

            gmic_use_argx;
            cimg_snprintf(argx,_argx.width(),"x = lerp(%g,%g,x/%d);",
                          xmin,xmax,(unsigned int)(resolution>1?resolution - 1:0));
            const CImg<char> n_formula = (CImg<char>::string(argx,false,true),
                                          CImg<char>::string(formula,true,true))>'x';
            boundary = 1U;
            try { // Determine vector dimension of specified formula
              typename CImg<T>::_cimg_math_parser mp(n_formula.data() + (*n_formula=='>' || *n_formula=='<' ||
                                                                         *n_formula=='*' || *n_formula==':'),
                                                     "graph",CImg<T>::const_empty(),0,&images);
              boundary = std::max(1U,mp.result_dim);
            } catch (...) { is_cond = false; }
            CImg<T> values((int)resolution,1,1,boundary,0);
            values.fill(n_formula,false,true,&images);

            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              g_img.assign(img.spectrum(),1,1,1,(T)0).fill_from_values(color,true);
              gmic_apply(gmic_draw_graph(values,g_img.data(),opacity,
                                         plot_type,vertex_type,ymin,ymax,pattern),true);
            }
          } else if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                                   gmic_use_indices,&sep,&end)==2 && sep==']') ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u%c",
                                  indices,&plot_type,&end)==2 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u%c",
                                  indices,&plot_type,&vertex_type,&end)==3 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%lf,%lf%c",
                                  indices,&plot_type,&vertex_type,&ymin,&ymax,&end)==5 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%lf,%lf,%f%c",
                                  indices,&plot_type,&vertex_type,&ymin,&ymax,&opacity,&end)==6||
                      (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%lf,%lf,%f,0%c%x%c",
                                   indices,&plot_type,&vertex_type,&ymin,&ymax,&opacity,&sep1,
                                   &pattern,&end)==8 &&
                       sep1=='x') ||
                      (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%lf,%lf,%f,"
                                   "%4095[0-9.eEinfa,+-]%c",
                                   indices,&plot_type,&vertex_type,&ymin,&ymax,&opacity,
                                   gmic_use_color,&end)==7 &&
                       (bool)(pattern=~0U)) ||
                      (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%lf,%lf,"
                                   "%f,0%c%x,%4095[0-9.eEinfa,+-]%c",
                                   indices,&plot_type,&vertex_type,&ymin,&ymax,
                                   &opacity,&sep1,&pattern,&(*color=0),&end)==9 &&
                       sep1=='x')) &&
                     (ind=selection2cimg(indices,images.size(),images_names,"graph")).height()==1 &&
                     plot_type<=3 && vertex_type<=7) {
            if (!plot_type && !vertex_type) plot_type = 1;
            print(images,0,"Draw graph of dataset [%u] on image%s, with %s contours, %s vertices, "
                  "y-range = (%g,%g), opacity %g, pattern 0x%x and color (%s).",
                  *ind,
                  gmic_selection.data(),
                  plot_type==0?"no":plot_type==1?"linear":plot_type==2?"spline":"bar",
                  vertex_type==0?"no":vertex_type==1?"dot":vertex_type==2?"straight cross":
                  vertex_type==3?"diagonal cross":vertex_type==4?"filled circle":
                  vertex_type==5?"outlined circle":vertex_type==6?"square":"diamond",
                  ymin,ymax,
                  opacity,pattern,
                  *color?color:"default");
            const CImg<T> values = gmic_image_arg(*ind);
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              g_img.assign(img.spectrum(),1,1,1,(T)0).fill_from_values(color,true);
              gmic_apply(draw_graph(values,g_img.data(),opacity,plot_type,vertex_type,ymin,ymax,pattern),true);
            }
          } else arg_error("graph");
          g_img.assign();
          is_change = true;
          ++position;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'h...'
        //-----------------------------
      gmic_commands_h :

        // Histogram.
        if (!std::strcmp("histogram",command)) {
          gmic_substitute_args(false);
          value = value0 = value1 = 0;
          sep = sep0 = sep1 = 0;
          is_cond = false; // no_min_max?
          if (((cimg_sscanf(argument,"%lf%c",
                            &value,&end)==1 && (is_cond = true)) ||
               ((cimg_sscanf(argument,"%lf%c%c",
                             &value,&sep,&end)==2 && sep=='%') && (is_cond = true)) ||
               cimg_sscanf(argument,"%lf,%lf,%lf%c",
                           &value,&value0,&value1,&end)==3 ||
               (cimg_sscanf(argument,"%lf%c,%lf,%lf%c",
                            &value,&sep,&value0,&value1,&end)==4 && sep=='%') ||
               (cimg_sscanf(argument,"%lf,%lf%c,%lf%c",
                            &value,&value0,&sep0,&value1,&end)==4 && sep0=='%') ||
               (cimg_sscanf(argument,"%lf%c,%lf%c,%lf%c",
                            &value,&sep,&value0,&sep0,&value1,&end)==5 && sep=='%' &&
                sep0=='%') ||
               (cimg_sscanf(argument,"%lf,%lf,%lf%c%c",
                            &value,&value0,&value1,&sep1,&end)==4 && sep1=='%') ||
               (cimg_sscanf(argument,"%lf%c,%lf,%lf%c%c",
                            &value,&sep,&value0,&value1,&sep1,&end)==5 && sep=='%' &&
                sep1=='%') ||
               (cimg_sscanf(argument,"%lf,%lf%c,%lf%c%c",
                            &value,&value0,&sep0,&value1,&sep1,&end)==5 && sep0=='%' &&
                sep1=='%') ||
               (cimg_sscanf(argument,"%lf%c,%lf%c,%lf%c%c",
                            &value,&sep,&value0,&sep0,&value1,&sep1,&end)==6 && sep=='%' &&
                sep0=='%' && sep1=='%')) &&
              value>=0.5) {
            value = cimg::round(value);
            if (is_cond) { value0 = 0; value1 = 100; sep0 = sep1 = '%'; }
            print(images,0,"Compute histogram of image%s, using %g%s level%s in range [%g%s,%g%s].",
                  gmic_selection.data(),
                  value,sep=='%'?"%":"",
                  value>1?"s":"",
                  value0,sep0=='%'?"%":"",
                  value1,sep1=='%'?"%":"");
            cimg_forY(selection,l) {
              CImg<T> &img = gmic_check(images[selection[l]]);
              nvalue0 = value0; nvalue1 = value1;
              vmin = vmax = 0;
              if (sep0=='%' || sep1=='%') {
                if (img) vmax = (double)img.max_min(vmin);
                if (sep0=='%') nvalue0 = vmin + (vmax - vmin)*value0/100;
                if (sep1=='%') nvalue1 = vmin + (vmax - vmin)*value1/100;
              }
              const unsigned int
                _nb_levels = std::max(1U,
                                      (unsigned int)cimg::round(sep=='%'?
                                                                value*(1 + nvalue1 - nvalue0)/100:
                                                                value));
              gmic_apply(histogram(_nb_levels,(T)nvalue0,(T)nvalue1),false);
            }
          } else arg_error("histogram");
          is_change = true;
          ++position;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'i...'
        //-----------------------------
      gmic_commands_i :
        if (is_command_input || // Redirect 'input'
            (command[1]=='f' && !command[2]) || // Redirect 'if'
            (command[1]=='f' && command[2]=='f' && command[3]=='t' && !command[4])) // Redirect 'ifft'
          goto gmic_commands_others;

        // Draw image.
        if (!std::strcmp("image",command)) {
          gmic_substitute_args(true);
          name.assign(256);
          double x = 0, y = 0, z = 0, c = 0;
          float max_opacity_mask = 1;
          *indices = *name = *argx = *argy = *argz = *argc = sep = sepx = sepy = sepz = sepc = 0;
          ind0.assign();
          opacity = 1;
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                            gmic_use_indices,&sep,&end)==2 && sep==']') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%~+-]%c",
                           indices,gmic_use_argx,&end)==2 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%~+-],%255[0-9.eE%~+-]%c",
                           indices,argx,gmic_use_argy,&end)==3 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%~+-],%255[0-9.eE%~+-],"
                           "%255[0-9.eE%~+-]%c",
                           indices,argx,argy,gmic_use_argz,&end)==4 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%~+-],%255[0-9.eE%~+-],"
                           "%255[0-9.eE%~+-],%255[0-9.eE%~+-]%c",
                           indices,argx,argy,argz,gmic_use_argc,&end)==5 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%~+-],%255[0-9.eE%~+-],"
                           "%255[0-9.eE%~+-],%255[0-9.eE%~+-],%f%c",
                           indices,argx,argy,argz,argc,&opacity,&end)==6 ||
               (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%~+-],%255[0-9.eE%~+-],"
                            "%255[0-9.eE%~+-],%255[0-9.eE%~+-],%f,[%255[a-zA-Z0-9_.%+-]%c%c",
                            indices,argx,argy,argz,argc,&opacity,name.data(),&sep,&end)==8 &&
                sep==']') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%~+-],%255[0-9.eE%~+-],"
                           "%255[0-9.eE%~+-],%255[0-9.eE%~+-],%f,[%255[a-zA-Z0-9_.%+-]],%f%c",
                           indices,argx,argy,argz,argc,&opacity,name.data(),
                           &max_opacity_mask,&end)==8) &&
              (ind=selection2cimg(indices,images.size(),images_names,"image")).height()==1 &&
              (!*name ||
               (ind0=selection2cimg(name,images.size(),images_names,"image")).height()==1) &&
              (!*argx ||
               cimg_sscanf(argx,"%lf%c",&x,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&x,&sepx,&end)==2 && (sepx=='%' || sepx=='~'))) &&
              (!*argy ||
               cimg_sscanf(argy,"%lf%c",&y,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&y,&sepy,&end)==2 && (sepy=='%' || sepy=='~'))) &&
              (!*argz ||
               cimg_sscanf(argz,"%lf%c",&z,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&z,&sepz,&end)==2 && (sepz=='%' || sepz=='~'))) &&
              (!*argc ||
               cimg_sscanf(argc,"%lf%c",&c,&end)==1 ||
               (cimg_sscanf(argc,"%lf%c%c",&c,&sepc,&end)==2 && (sepc=='%' || sepc=='~')))) {
            const CImg<T> sprite = gmic_image_arg(*ind);
            CImg<T> mask;
            if (ind0) {
              mask = gmic_image_arg(*ind0);
              print(images,0,"Draw image [%u] at (%g%s,%g%s,%g%s,%g%s) on image%s, "
                    "with opacity %g and mask [%u].",
                    *ind,
                    x,sepx=='%'?"%":sepx=='~'?"~":"",
                    y,sepy=='%'?"%":sepy=='~'?"~":"",
                    z,sepz=='%'?"%":sepz=='~'?"~":"",
                    c,sepc=='%'?"%":sepc=='~'?"~":"",
                    gmic_selection.data(),
                    opacity,
                    *ind0);
            } else print(images,0,"Draw image [%u] at (%g%s,%g%s,%g%s,%g%s) on image%s, "
                         "with opacity %g.",
                         *ind,
                         x,sepx=='%'?"%":sepx=='~'?"~":"",
                         y,sepy=='%'?"%":sepy=='~'?"~":"",
                         z,sepz=='%'?"%":sepz=='~'?"~":"",
                         c,sepc=='%'?"%":sepc=='~'?"~":"",
                         gmic_selection.data(),
                         opacity);
            cimg_forY(selection,l)
              if (ind0) {
                gmic_apply(gmic_draw_image((float)x,(float)y,(float)z,(float)c,sepx,sepy,sepz,sepc,
                                           sprite,mask,opacity,max_opacity_mask),true);
              }
              else {
                gmic_apply(gmic_draw_image((float)x,(float)y,(float)z,(float)c,
                                           sepx,sepy,sepz,sepc,sprite,opacity),true);
              }
          } else arg_error("image");
          is_change = true;
          ++position;
          continue;
        }

        // Index image with a LUT.
        if (!std::strcmp("index",command)) {
          gmic_substitute_args(true);
          unsigned int map_indexes = 0;
          float dithering = 0;
          sep = 0;
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                            gmic_use_indices,&sep,&end)==2 && sep==']') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%f%c",
                           indices,&dithering,&end)==2 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%f,%u%c",
                           indices,&dithering,&map_indexes,&end)==3) &&
              (ind=selection2cimg(indices,images.size(),images_names,"index")).height()==1) {
            const float ndithering = dithering<0?0:dithering>1?1:dithering;
            print(images,0,"Index values in image%s by LUT [%u], with dithering level %g%s.",
                  gmic_selection.data(),
                  *ind,
                  ndithering,
                  map_indexes?" and index mapping":"");
            const CImg<T> palette = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply(index(palette,ndithering,(bool)map_indexes),false);
            is_change = true;
            ++position;
            continue;
          }
          // When command 'index' is invoked with different arguments, custom version in stdlib
          // is used rather than the built-in version.
          is_builtin_command = false;
          goto gmic_commands_others;
        }

        // Matrix inverse (or pseudoinverse).
        if (!std::strcmp("invert",command)) {
          gmic_substitute_args(false);
          pattern = 0; value = 0;
          if ((cimg_sscanf(argument,"%u%c",
                           &pattern,&end)==1 ||
               cimg_sscanf(argument,"%u,%lf%c",
                           &pattern,&value,&end)==2) &&
              pattern<=1 && value>=0) ++position;
          else pattern = 0;
          print(images,0,"Invert matrix image%s, using %s solver and lambda %g.",
                gmic_selection.data(),pattern?"LU":"SVD",value);
          cimg_forY(selection,l) gmic_apply(invert((bool)pattern,(float)value),false);
          is_change = true;
          continue;
        }

        // Extract 3D isoline.
        if (!std::strcmp("isoline3d",command)) {
          gmic_substitute_args(false);
          double x0 = -3, y0 = -3, x1 = 3, y1 = 3, dx = 256, dy = 256;
          sep = sepx = sepy = *formula = 0;
          value = 0;
          if (cimg_sscanf(argument,"%lf%c",
                          &value,&end)==1 ||
              cimg_sscanf(argument,"%lf%c%c",
                          &value,&sep,&end)==2) {
            print(images,0,"Extract 3D isolines from image%s, using isovalue %g%s.",
                  gmic_selection.data(),
                  value,sep=='%'?"%":"");
            cimg_forY(selection,l) {
              uind = selection[l];
              CImg<T>& img = gmic_check(images[uind]);
              if (img) {
                vertices.assign();
                primitives.assign();
                g_list_uc.assign();
                g_img_uc.assign(3,img.spectrum(),1,1,220).noise(35,1);
                if (img.spectrum()==1) g_img_uc(0) = g_img_uc(1) = g_img_uc(2) = 200;
                else {
                  g_img_uc(0,0) = 255; g_img_uc(1,0) = g_img_uc(2,0) = 30;
                  g_img_uc(0,1) = g_img_uc(2,1) = 30; g_img_uc(1,1) = 255;
                  if (img.spectrum()>=3) { g_img_uc(0,2) = g_img_uc(1,2) = 30; g_img_uc(2,2) = 255; }
                }
                cimg_forC(img,k) {
                  const CImg<T> channel = img.get_shared_channel(k);
                  nvalue = value;
                  if (sep=='%') {
                    vmax = (double)channel.max_min(vmin);
                    nvalue = vmin + (vmax - vmin)*value/100;
                  }
                  CImgList<unsigned int> prims;
                  const CImg<float> pts = img.get_shared_channel(k).get_isoline3d(prims,(float)nvalue);
                  vertices.append_object3d(primitives,pts,prims);
                  g_list_uc.insert(prims.size(),CImg<unsigned char>::vector(g_img_uc(0,k),
                                                                            g_img_uc(1,k),
                                                                            g_img_uc(2,k)));
                }
                if (!vertices)
                  warn(images,0,false,
                       "Command 'isoline3d': Isovalue %g%s not found in image [%u].",
                       value,sep=='%'?"%":"",uind);
                vertices.object3dtoCImg3d(primitives,g_list_uc,false);
                gmic_apply(replace(vertices),false);
                primitives.assign();
                g_list_uc.assign();
                g_img_uc.assign();
              } else gmic_apply(replace(img),false);
            }
          } else if ((cimg_sscanf(argument,"'%4095[^']',%lf%c",
                                  gmic_use_formula,&value,&end)==2 ||
                      cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf%c",
                                  formula,&value,&x0,&y0,&x1,&y1,&end)==6 ||
                      cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf,%lf%c",
                                  formula,&value,&x0,&y0,&x1,&y1,&dx,&dy,&end)==8 ||
                      (cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf%c,%lf%c",
                                   formula,&value,&x0,&y0,&x1,&y1,&dx,&sepx,&dy,&end)==9 &&
                       sepx=='%') ||
                      (cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf,%lf%c%c",
                                   formula,&value,&x0,&y0,&x1,&y1,&dx,&dy,&sepy,&end)==9 &&
                       sepy=='%') ||
                      (cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf%c,%lf%c%c",
                                   formula,&value,&x0,&y0,&x1,&y1,&dx,&sepx,&dy,&sepy,&end)==10&&
                       sepx=='%' && sepy=='%')) &&
                     dx>0 && dy>0) {
            dx = cimg::round(dx);
            dy = cimg::round(dy);
            strreplace_fw(formula);
            print(images,0,"Extract 3D isoline %g from formula '%s', in range (%g,%g)-(%g,%g) "
                  "with size %g%sx%g%s.",
                  value,
                  formula,
                  x0,y0,
                  x1,y1,
                  dx,sepx=='%'?"%":"",
                  dy,sepy=='%'?"%":"");
            if (sepx=='%') dx = -dx;
            if (sepy=='%') dy = -dy;
            CImg<T>::isoline3d(primitives,(const char*)formula,(float)value,
                               (float)x0,(float)y0,(float)x1,(float)y1,(int)dx,(int)dy).move_to(vertices);
            vertices.object3dtoCImg3d(primitives,false).move_to(images);
            primitives.assign();
            gmic_use_title;
            cimg_snprintf(title,_title.width(),"[3D isoline %g of '%s']",value,formula);
            CImg<char>::string(title).move_to(images_names);
          } else arg_error("isoline3d");
          is_change = true;
          ++position;
          continue;
        }

        // Extract 3D isosurface.
        if (!std::strcmp("isosurface3d",command)) {
          gmic_substitute_args(false);
          double x0 = -3, y0 = -3, z0 = -3, x1 = 3, y1 = 3, z1 = 3,
            dx = 32, dy = 32, dz = 32;
          sep = sepx = sepy = sepz = *formula = 0;
          value = 0;
          if (cimg_sscanf(argument,"%lf%c",
                          &value,&end)==1 ||
              cimg_sscanf(argument,"%lf%c%c",
                          &value,&sep,&end)==2) {
            print(images,0,"Extract 3D isosurface from image%s, using isovalue %g%s.",
                  gmic_selection.data(),
                  value,sep=='%'?"%":"");
            cimg_forY(selection,l) {
              uind = selection[l];
              CImg<T>& img = gmic_check(images[uind]);
              if (img) {
                vertices.assign();
                primitives.assign();
                g_list_uc.assign();
                g_img_uc.assign(3,img.spectrum(),1,1,220).noise(35,1);
                if (img.spectrum()==1) g_img_uc(0) = g_img_uc(1) = g_img_uc(2) = 200;
                else {
                  g_img_uc(0,0) = 255; g_img_uc(1,0) = g_img_uc(2,0) = 30;
                  g_img_uc(0,1) = g_img_uc(2,1) = 30; g_img_uc(1,1) = 255;
                  if (img.spectrum()>=3) { g_img_uc(0,2) = g_img_uc(1,2) = 30; g_img_uc(2,2) = 255; }
                }
                cimg_forC(img,k) {
                  const CImg<T> channel = img.get_shared_channel(k);
                  nvalue = value;
                  if (sep=='%') {
                    vmax = (double)channel.max_min(vmin);
                    nvalue = vmin + (vmax - vmin)*value/100;
                  }
                  CImgList<unsigned int> prims;
                  const CImg<float> pts = channel.get_isosurface3d(prims,(float)nvalue);
                  vertices.append_object3d(primitives,pts,prims);
                  g_list_uc.insert(prims.size(),CImg<unsigned char>::vector(g_img_uc(0,k),
                                                                            g_img_uc(1,k),
                                                                            g_img_uc(2,k)));
                }
                if (!vertices) {
                  if (img.depth()>1)
                    warn(images,0,false,
                         "Command 'isosurface3d': Isovalue %g%s not found in image [%u].",
                         value,sep=='%'?"%":"",uind);
                  else
                    warn(images,0,false,
                         "Command 'isosurface3d': Image [%u] has a single slice, "
                         "isovalue %g%s not found.",
                         uind,value,sep=='%'?"%":"");
                }
                vertices.object3dtoCImg3d(primitives,g_list_uc,false);
                gmic_apply(replace(vertices),false);
                primitives.assign(); g_list_uc.assign(); g_img_uc.assign();
              } else gmic_apply(replace(img),false);
            }
          } else if ((cimg_sscanf(argument,"'%4095[^']',%lf%c",
                                  gmic_use_formula,&value,&end)==2 ||
                      cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf,%lf%c",
                                  formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,&end)==8 ||
                      cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c",
                                  formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,&dx,&dy,&dz,&end)==11 ||
                      (cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c,%lf,%lf%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,
                                   &dx,&sepx,&dy,&dz,&end)==12 &&
                       sepx=='%') ||
                      (cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c,%lf%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,
                                   &dx,&dy,&sepy,&dz,&end)==12 &&
                       sepy=='%') ||
                      (cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,
                                   &dx,&dy,&dz,&sepz,&end)==12 &&
                       sepz=='%') ||
                      (cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c,%lf%c,%lf%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,
                                   &dx,&sepx,&dy,&sepy,&dz,&end)==13 &&
                       sepx=='%' && sepy=='%') ||
                      (cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c,%lf,%lf%c%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,
                                   &dx,&sepx,&dy,&dz,&sepz,&end)==13 &&
                       sepx=='%' && sepz=='%') ||
                      (cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c,%lf%c%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,
                                   &dx,&dy,&sepy,&dz,&sepz,&end)==13 &&
                       sepy=='%' && sepz=='%') ||
                      (cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c,%lf%c,%lf%c%c",
                                   formula,&value,&x0,&y0,&z0,&x1,&y1,&z1,
                                   &dx,&sepx,&dy,&sepy,&dz,&sepz,&end)==14 &&
                       sepx=='%' && sepy=='%' && sepz=='%')) &&
                     dx>0 && dy>0 && dz>0) {
            dx = cimg::round(dx);
            dy = cimg::round(dy);
            dz = cimg::round(dz);
            strreplace_fw(formula);
            print(images,0,"Extract 3D isosurface %g from formula '%s', "
                  "in range (%g,%g,%g)-(%g,%g,%g) with size %g%sx%g%sx%g%s.",
                  value,
                  formula,
                  x0,y0,z0,
                  x1,y1,z1,
                  dx,sepx=='%'?"%":"",
                  dy,sepy=='%'?"%":"",
                  dz,sepz=='%'?"%":"");
            if (sepx=='%') dx = -dx;
            if (sepy=='%') dy = -dy;
            if (sepz=='%') dz = -dz;
            CImg<T>::isosurface3d(primitives,(const char*)formula,(float)value,
                                  (float)x0,(float)y0,(float)z0,(float)x1,(float)y1,(float)z1,
                                  (int)dx,(int)dy,(int)dz).move_to(vertices);
            vertices.object3dtoCImg3d(primitives,false).move_to(images);
            primitives.assign();
            gmic_use_title;
            cimg_snprintf(title,_title.width(),"[3D isosurface %g of '%s']",value,formula);
            CImg<char>::string(title).move_to(images_names);
          } else arg_error("isosurface3d");
          is_change = true;
          ++position;
          continue;
        }

        // Inpaint.
        if (!std::strcmp("inpaint",command)) {
          gmic_substitute_args(true);
          double patch_size = 11, lookup_size = 22, lookup_factor = 0.5, lookup_increment = 1,
            blend_size = 0, blend_threshold = 0, blend_decay = 0.05f, blend_scales = 10;
          unsigned int is_blend_outer = 1, method = 1;
          sep = *indices = 0;
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep,&end)==2 &&
                sep==']') ||
               (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%c%c",indices,&sep,&end)==2 &&
                sep=='0') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],0,%u%c",indices,&method,&end)==2) &&
              (ind=selection2cimg(indices,images.size(),images_names,"inpaint")).height()==1 &&
              method<=3) {
            print(images,0,"Inpaint image%s masked by image [%u], with %s algorithm.",
                  gmic_selection.data(),
                  *ind,
                  method==0?"low-connectivity average":method==1?"high-connectivity average":
                  method==2?"low-connectivity median":"high-connectivity median");
            const CImg<T> mask = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply(inpaint(mask,method),false);
          } else if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                                   gmic_use_indices,&sep,&end)==2 && sep==']') ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf%c",
                                  indices,&patch_size,&end)==2 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf%c",
                                  indices,&patch_size,&lookup_size,&end)==3 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf%c",
                                  indices,&patch_size,&lookup_size,&lookup_factor,&end)==4 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf,%lf%c",
                                  indices,&patch_size,&lookup_size,&lookup_factor,
                                  &lookup_increment,&end)==5 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf,%lf,%lf%c",
                                  indices,&patch_size,&lookup_size,&lookup_factor,
                                  &lookup_increment,&blend_size,&end)==6 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf,%lf,%lf,%lf%c",
                                  indices,&patch_size,&lookup_size,&lookup_factor,
                                  &lookup_increment,&blend_size,&blend_threshold,&end)==7 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf,%lf,%lf,%lf,%lf%c",
                                  indices,&patch_size,&lookup_size,&lookup_factor,
                                  &lookup_increment,&blend_size,&blend_threshold,&blend_decay,
                                  &end)==8 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c",
                                  indices,&patch_size,&lookup_size,&lookup_factor,
                                  &lookup_increment,&blend_size,&blend_threshold,&blend_decay,
                                  &blend_scales,&end)==9 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%u%c",
                                  indices,&patch_size,&lookup_size,&lookup_factor,
                                  &lookup_increment,&blend_size,&blend_threshold,&blend_decay,
                                  &blend_scales,&is_blend_outer,&end)==10) &&
                     (ind=selection2cimg(indices,images.size(),images_names,"inpaint")).height()==1 &&
                     patch_size>=0.5 && lookup_size>=0.5 && lookup_factor>=0 &&
                     blend_size>=0 && blend_threshold>=0 && blend_threshold<=1 &&
                     blend_decay>=0 && blend_scales>=0.5 && is_blend_outer<=1) {
            const CImg<T> mask = gmic_image_arg(*ind);
            patch_size = cimg::round(patch_size);
            lookup_size = cimg::round(lookup_size);
            lookup_increment = cimg::round(lookup_increment);
            blend_size = cimg::round(blend_size);
            blend_scales = cimg::round(blend_scales);
            print(images,0,"Inpaint image%s masked by image [%d], with patch size %g, "
                  "lookup size %g, lookup factor %g, lookup_increment %g, blend size %g, "
                  "blend threshold %g, blend decay %g, %g blend scale%s and outer blending %s.",
                  gmic_selection.data(),*ind,
                  patch_size,lookup_size,lookup_factor,lookup_increment,
                  blend_size,blend_threshold,blend_decay,blend_scales,blend_scales!=1?"s":"",
                  is_blend_outer?"enabled":"disabled");
            cimg_forY(selection,l)
              gmic_apply(inpaint_patch(mask,
                                       (unsigned int)patch_size,(unsigned int)lookup_size,
                                       (float)lookup_factor,
                                       (int)lookup_increment,
                                       (unsigned int)blend_size,(float)blend_threshold,(float)blend_decay,
                                       (unsigned int)blend_scales,(bool)is_blend_outer),false);
          } else arg_error("inpaint");
          is_change = true;
          ++position;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'k...'
        //-----------------------------
      gmic_commands_k :

        // Keep images.
        if (!std::strcmp("keep",command)) {
          print(images,0,"Keep image%s",
                gmic_selection.data());
          g_list.assign(selection.height());
          g_list_c.assign(selection.height());
          if (is_get) {
            cimg_forY(selection,l) {
              uind = selection[l];
              g_list[l].assign(images[uind]);
              images_names[uind].get_copymark().move_to(g_list_c[l]);
            }
            g_list.move_to(images,~0U);
            g_list_c.move_to(images_names,~0U);
          } else {
            cimg_forY(selection,l) {
              uind = selection[l];
              g_list[l].swap(images[uind]);
              g_list_c[l].swap(images_names[uind]);
            }
            g_list.swap(images);
            g_list_c.swap(images_names);
          }
          if (is_verbose) {
            cimg::mutex(29);
            std::fprintf(cimg::output()," (%u image%s left).",
                         images.size(),images.size()==1?"":"s");
            std::fflush(cimg::output());
            cimg::mutex(29,0);
          }
          g_list.assign();
          g_list_c.assign();
          is_change = true;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'l...'
        //-----------------------------
      gmic_commands_l :

        // Start local environment.
        if (!std::strcmp("local",command)) {
          if (is_debug_info && debug_line!=~0U) {
            gmic_use_argx;
            cimg_snprintf(argx,_argx.width(),"*local#%u",debug_line);
            CImg<char>::string(argx).move_to(callstack);
          } else CImg<char>::string("*local").move_to(callstack);
          if (is_very_verbose)
            print(images,0,"Start 'local...done' block, with image%s.",
                  gmic_selection.data());
          g_list.assign(selection.height());
          g_list_c.assign(selection.height());

          if (is_get) cimg_forY(selection,l) {
              uind = selection[l];
              g_list[l].assign(images[uind]);
              g_list_c[l] = images_names[uind];
            } else {
            cimg::mutex(27);
            cimg_forY(selection,l) {
              uind = selection[l];
              if (images[uind].is_shared())
                g_list[l].assign(images[uind],false);
              else {
                if ((images[uind].width() || images[uind].height()) && !images[uind]._spectrum) {
                  selection2string(selection,images_names,1,name);
                  error(true,images,0,0,
                        "Command 'local': Invalid selection%s "
                        "(image [%u] is already used in another thread).",
                        name.data() + (*name=='s'?1:0),uind);
                }
                g_list[l].swap(images[uind]);
                // Small hack to be able to track images of the selection passed to the new environment.
                std::memcpy(&images[uind]._width,&g_list[l]._data,sizeof(void*));
                images[uind]._spectrum = 0;
              }
              g_list_c[l] = images_names[uind]; // Make a copy to be still able to recognize 'pass[label]'
            }
            cimg::mutex(27,0);
          }

          is_lbrace_command = true;
          const int o_verbosity = verbosity;
          gmic_exception exception;
          try {
            if (next_debug_line!=~0U) { debug_line = next_debug_line; next_debug_line = ~0U; }
            if (next_debug_filename!=~0U) { debug_filename = next_debug_filename; next_debug_filename = ~0U; }
            _run(commands_line,position,g_list,g_list_c,images,images_names,variables_sizes,is_noarg,0,
                 command_selection,false);
          } catch (gmic_exception &e) {
            check_elif = false;
            int nb_levels = 0;
            for (nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
              it = commands_line[position];
              if (*it==1)
                is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
              else {
                _is_get = *it=='+';
                if (*it==1)
                  is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
                else {
                  it+=(_is_get || *it=='-');
                  gmic_if_flr ++nb_levels; gmic_elif_flr --nb_levels;
                  else if (!_is_get && nb_levels==1 && !std::strcmp("onfail",it)) break;
                }
              }
            }
            if (nb_levels==1 && position<commands_line.size()) { // Onfail block found
              verbosity = o_verbosity; // Restore verbosity
              if (is_very_verbose) print(images,0,"Reach 'onfail' block.");
              try {
                _run(commands_line,++position,g_list,g_list_c,
                     parent_images,parent_images_names,variables_sizes,is_noarg,0,0,false);
              } catch (gmic_exception &e2) {
                cimg::swap(exception._command,e2._command);
                cimg::swap(exception._message,e2._message);
              }
            } else {
              cimg::swap(exception._command,e._command);
              cimg::swap(exception._message,e._message);
            }
          }

          cimg::mutex(27);
          if (is_get) {
            g_list.move_to(images,~0U);
            cimglist_for(g_list_c,i) g_list_c[i].get_copymark().move_to(g_list_c[i]);
            g_list_c.move_to(images_names,~0U);
          } else {
            const unsigned int nb = std::min((unsigned int)selection.height(),g_list.size());
            if (nb>0) {
              for (unsigned int i = 0; i<nb; ++i) {
                uind = selection[i];
                if (images[uind].is_shared()) {
                  images[uind] = g_list[i];
                  g_list[i].assign();
                } else images[uind].swap(g_list[i]);
                images_names[uind].swap(g_list_c[i]);
              }
              g_list.remove(0,nb - 1);
              g_list_c.remove(0,nb - 1);
            }
            if (nb<(unsigned int)selection.height())
              remove_images(images,images_names,selection,nb,selection.height() - 1);
            else if (g_list) {
              const unsigned int uind0 = selection?selection.back() + 1:images.size();
              images.insert(g_list,uind0);
              g_list_c.move_to(images_names,uind0);
            }
          }
          g_list.assign();
          g_list_c.assign();
          if (!exception._message) callstack.remove();
          cimg::mutex(27,0);
          if (exception._message) throw exception;
          continue;
        }

        // Less or equal.
        gmic_arithmetic_command("le",
                                operator_le,
                                "Compute boolean 'less or equal than' between image%s and %g%s",
                                gmic_selection.data(),value,ssep,T,
                                operator_le,
                                "Compute boolean 'less or equal than' between image%s and image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_le,
                                "Compute boolean 'less or equal than' between image%s and "
                                "expression %s'",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute boolean 'less or equal than' between image%s");

        // Less than.
        gmic_arithmetic_command("lt",
                                operator_lt,
                                "Compute boolean 'less than' between image%s and %g%s",
                                gmic_selection.data(),value,ssep,T,
                                operator_lt,
                                "Compute boolean 'less than' between image%s and image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_lt,
                                "Compute boolean 'less than' between image%s and expression %s'",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute boolean 'less than' between image%s");

        // Logarithm, base-e.
        gmic_simple_command("log",log,"Compute pointwise base-e logarithm of image%s.");

        // Logarithm, base-2.
        gmic_simple_command("log2",log2,"Compute pointwise base-2 logarithm of image%s.");

        // Logarithm, base-10.
        gmic_simple_command("log10",log10,"Compute pointwise base-10 logarithm of image%s.");

        // Draw line.
        if (!std::strcmp("line",command)) {
          gmic_substitute_args(false);
          *argx = *argy = *argz = *argc = *color = 0;
          double x0 = 0, y0 = 0, x1 = 0, y1 = 0;
          char sepx0 = 0, sepy0 = 0, sepx1 = 0, sepy1 = 0;
          sep1 = 0;
          pattern = ~0U; opacity = 1;
          if ((cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-]%c",
                           gmic_use_argx,gmic_use_argy,gmic_use_argz,gmic_use_argc,&end)==4 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-],%f%c",
                           argx,argy,argz,argc,&opacity,&end)==5 ||
               (cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                            "%255[0-9.eE%+-],%f,0%c%x%c",
                            argx,argy,argz,argc,&opacity,&sep1,&pattern,&end)==7 &&
                sep1=='x') ||
               (cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                            "%255[0-9.eE%+-],%f,%4095[0-9.eEinfa,+-]%c",
                            argx,argy,argz,argc,&opacity,gmic_use_color,&end)==6 && (bool)(pattern=~0U)) ||
               (cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                            "%255[0-9.eE%+-],%f,0%c%x,%4095[0-9.eEinfa,+-]%c",
                            argx,argy,argz,argc,&opacity,&sep1,
                            &pattern,&(*color=0),&end)==8 && sep1=='x')) &&
              (cimg_sscanf(argx,"%lf%c",&x0,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&x0,&sepx0,&end)==2 && sepx0=='%')) &&
              (cimg_sscanf(argy,"%lf%c",&y0,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&y0,&sepy0,&end)==2 && sepy0=='%')) &&
              (cimg_sscanf(argz,"%lf%c",&x1,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&x1,&sepx1,&end)==2 && sepx1=='%')) &&
              (cimg_sscanf(argc,"%lf%c",&y1,&end)==1 ||
               (cimg_sscanf(argc,"%lf%c%c",&y1,&sepy1,&end)==2 && sepy1=='%'))) {
            print(images,0,"Draw line (%g%s,%g%s) - (%g%s,%g%s) on image%s, with opacity %g, "
                  "pattern 0x%x and color (%s).",
                  x0,sepx0=='%'?"%":"",
                  y0,sepy0=='%'?"%":"",
                  x1,sepx1=='%'?"%":"",
                  y1,sepy1=='%'?"%":"",
                  gmic_selection.data(),
                  opacity,pattern,
                  *color?color:"default");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              g_img.assign(img.spectrum(),1,1,1,(T)0).fill_from_values(color,true);
              const int
                nx0 = (int)cimg::round(sepx0=='%'?x0*(img.width() - 1)/100:x0),
                ny0 = (int)cimg::round(sepy0=='%'?y0*(img.height() - 1)/100:y0),
                nx1 = (int)cimg::round(sepx1=='%'?x1*(img.width() - 1)/100:x1),
                ny1 = (int)cimg::round(sepy1=='%'?y1*(img.height() - 1)/100:y1);
              gmic_apply(draw_line(nx0,ny0,nx1,ny1,g_img.data(),opacity,pattern),true);
            }
          } else arg_error("line");
          g_img.assign();
          is_change = true;
          ++position;
          continue;
        }

        // Label connected components.
        if (!std::strcmp("label",command)) {
          gmic_substitute_args(false);
          float tolerance = 0;
          unsigned int is_L2_norm = 1;
          is_high_connectivity = 0;
          if ((cimg_sscanf(argument,"%f%c",&tolerance,&end)==1 ||
               cimg_sscanf(argument,"%f,%u%c",&tolerance,&is_high_connectivity,&end)==2 ||
               cimg_sscanf(argument,"%f,%u,%u%c",&tolerance,&is_high_connectivity,&is_L2_norm,&end)==3) &&
              tolerance>=0 && is_high_connectivity<=1 && is_L2_norm<=1) ++position;
          else { tolerance = 0; is_high_connectivity = 0; is_L2_norm = 1; }
          print(images,0,
                "Label connected components on image%s, with tolerance %g (L%d-norm) and "
                "%s connectivity.",
                gmic_selection.data(),tolerance,1 + is_L2_norm,is_high_connectivity?"high":"low");
          cimg_forY(selection,l) gmic_apply(label((bool)is_high_connectivity,tolerance,is_L2_norm),false);
          is_change = true;
          continue;
        }

        // Set 3D light position.
        if (!is_get && !std::strcmp("light3d",item)) {
          gmic_substitute_args(true);
          float lx = 0, ly = 0, lz = -5e8f;
          sep = *indices = 0;
          if (cimg_sscanf(argument,"%f,%f,%f%c",
                          &lx,&ly,&lz,&end)==3) {
            print(images,0,"Set 3D light position to (%g,%g,%g).",
                  lx,ly,lz);
            light3d_x = lx;
            light3d_y = ly;
            light3d_z = lz;
            ++position;
          } else if (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep,&end)==2 &&
                     sep==']' &&
                     (ind=selection2cimg(indices,images.size(),images_names,"light3d")).height()==1) {
            print(images,0,"Set 3D light texture from image [%u].",*ind);
            light3d.assign(images[*ind],false);
            ++position;
          } else {
            print(images,0,"Reset 3D light to default.");
            light3d.assign();
            light3d_x = light3d_y = 0; light3d_z = -5e8f;
          }
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'm...'
        //-----------------------------
      gmic_commands_m :
        if (command[1]=='u' && command[2]=='l' && command[3]=='3' && command[4]=='d' && !command[5]) // Redirect 'mul3d'
          goto gmic_commands_others;

        // Move images.
        if (!std::strcmp("move",command)) {
          gmic_substitute_args(false);
          double pos = 0;
          sep = 0;
          if (cimg_sscanf(argument,"%lf%c",&pos,&end)==1 ||
              (cimg_sscanf(argument,"%lf%c%c",&pos,&sep,&end)==2 && sep=='%')) {
            const int
              _iind0 = (int)cimg::round(sep=='%'?pos*images.size()/100:pos),
              iind0 = _iind0<0?_iind0 + (int)images.size():_iind0;
            if (iind0<0 || iind0>(int)images.size())
              error(true,images,0,0,
                    "Command 'move': Invalid position %d (not in range -%u...%u).",
                    _iind0,images.size(),images.size() - 1);
            print(images,0,"Move image%s to position %d.",
                  gmic_selection.data(),
                  iind0);
            CImgList<T> _images, nimages;
            CImgList<char> _images_names, nimages_names;
            if (is_get) {
              _images.insert(images.size());
              // Copy original list while preserving shared state of each item.
              cimglist_for(_images,l) _images[l].assign(images[l],images[l].is_shared());
              _images_names.assign(images_names);
            }
            nimages.insert(selection.height());
            cimg_forY(selection,l) {
              uind = selection[l];
              if (is_get) images[uind].move_to(nimages[l]);
              else images[uind].swap(nimages[l]);
              // Empty shared image as a special item to be removed later.
              images[uind]._is_shared = true;
              images_names[uind].move_to(nimages_names);
            }
            images.insert(nimages.size(),iind0);
            cimglist_for(nimages,l) nimages[l].swap(images[iind0 + l]);
            nimages_names.move_to(images_names,iind0);
            cimglist_for(images,l) if (!images[l] && images[l].is_shared()) { // Remove special items
              images.remove(l);
              images_names.remove(l--);
            }
            if (is_get) {
              cimglist_for(images,l) { // Replace shared items by non-shared ones for a get version
                if (images[l].is_shared()) {
                  CImg<T> tmp;
                  (images[l].move_to(tmp)).swap(images[l]);
                }
                images_names[l].get_copymark().move_to(images_names[l]);
              }
              images.insert(_images.size(),0);
              cimglist_for(_images,l) images[l].swap(_images[l]);
              _images_names.move_to(images_names,0);
            }
          } else arg_error("move");
          is_change = true;
          ++position;
          continue;
        }

        // Multiplication.
        gmic_arithmetic_command("mul",
                                operator*=,
                                "Multiply image%s by %g%s",
                                gmic_selection.data(),value,ssep,Tfloat,
                                mul,
                                "Multiply image%s by image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_muleq,
                                "Multiply image%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Multiply image%s");
        // Modulo.
        gmic_arithmetic_command("mod",
                                operator%=,
                                "Compute pointwise modulo of image%s by %g%s",
                                gmic_selection.data(),value,ssep,T,
                                operator%=,
                                "Compute pointwise modulo of image%s by image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_modeq,
                                "Compute pointwise modulo of image%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute sequential pointwise modulo of image%s");

        // Max.
        gmic_arithmetic_command("max",
                                max,
                                "Compute pointwise maximum between image%s and %g%s",
                                gmic_selection.data(),value,ssep,T,
                                max,
                                "Compute pointwise maximum between image%s and image [%d]",
                                gmic_selection.data(),ind[0],
                                max,
                                "Compute pointwise maximum between image%s and expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute pointwise maximum of all image%s together");

        // Maxabs.
        gmic_arithmetic_command("maxabs",
                                maxabs,
                                "Compute pointwise maxabs between image%s and %g%s",
                                gmic_selection.data(),value,ssep,T,
                                maxabs,
                                "Compute pointwise maxabs between image%s and image [%d]",
                                gmic_selection.data(),ind[0],
                                maxabs,
                                "Compute pointwise maxabs between image%s and expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute pointwise maxabs of all image%s together");

        // Min.
        gmic_arithmetic_command("min",
                                min,
                                "Compute pointwise minimum between image%s and %g%s",
                                gmic_selection.data(),value,ssep,T,
                                min,
                                "Compute pointwise minimum between image%s and image [%d]",
                                gmic_selection.data(),ind[0],
                                min,
                                "Compute pointwise minimum between image%s and expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute pointwise minimum of image%s");

        // Minabs.
        gmic_arithmetic_command("minabs",
                                minabs,
                                "Compute pointwise minabs between image%s and %g%s",
                                gmic_selection.data(),value,ssep,T,
                                minabs,
                                "Compute pointwise minabs between image%s and image [%d]",
                                gmic_selection.data(),ind[0],
                                minabs,
                                "Compute pointwise minabs between image%s and expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute pointwise minabs of image%s");

        // Matrix multiplication.
        gmic_arithmetic_command("mmul",
                                operator*=,
                                "Multiply matrix/vector%s by %g%s",
                                gmic_selection.data(),value,ssep,Tfloat,
                                operator*=,
                                "Multiply matrix/vector%s by matrix/vector image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_muleq,
                                "Multiply matrix/vector%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Multiply matrix/vector%s");

        // Map LUT.
        if (!std::strcmp("map",command)) {
          gmic_substitute_args(true);
          sep = *indices = 0;
          boundary = 0;
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep,&end)==2 && sep==']') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u%c",indices,&boundary,&end)==2) &&
              (ind=selection2cimg(indices,images.size(),images_names,"map")).height()==1 &&
              boundary<=3) {
            print(images,0,"Map LUT [%u] on image%s, with %s boundary conditions.",
                  *ind,
                  gmic_selection.data(),
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror");
            const CImg<T> palette = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply(map(palette,boundary),false);
            is_change = true;
            ++position;
            continue;
          }
          // If command 'map' is invoked with different arguments, custom version in stdlib
          // is used rather than the built-in version.
          is_builtin_command = false;
          goto gmic_commands_others;
        }

        // Mirror.
        if (!std::strcmp("mirror",command)) {
          gmic_substitute_args(false);
          bool is_valid_argument = *argument!=0;
          if (is_valid_argument) for (const char *s = argument; *s; ++s) {
              const char _s = *s;
              if (_s!='x' && _s!='y' && _s!='z' && _s!='c') { is_valid_argument = false; break; }
            }
          if (is_valid_argument) {
            print(images,0,"Mirror image%s along the '%s'-ax%cs.",
                  gmic_selection.data(),
                  gmic_argument_text_printed(),
                  std::strlen(argument)>1?'e':'i');
            cimg_forY(selection,l) gmic_apply(mirror(argument),false);
          } else arg_error("mirror");
          is_change = true;
          ++position;
          continue;
        }

        // Median filter.
        if (!std::strcmp("median",command)) {
          gmic_substitute_args(false);
          double fsiz = 3;
          float threshold = 0;
          if ((cimg_sscanf(argument,"%lf%c",
                           &fsiz,&end)==1 ||
               cimg_sscanf(argument,"%lf,%f%c",
                           &fsiz,&threshold,&end)==2) &&
              fsiz>=0 && threshold>=0) {
            fsiz = cimg::round(fsiz);
            if (threshold)
              print(images,0,"Apply median filter of size %g with threshold %g, on image%s.",
                    fsiz,threshold,
                    gmic_selection.data());
            else
              print(images,0,"Apply median filter of size %g, on image%s.",
                    fsiz,
                    gmic_selection.data());
            cimg_forY(selection,l) gmic_apply(blur_median((unsigned int)fsiz,threshold),false);
          } else arg_error("median");
          is_change = true;
          ++position;
          continue;
        }

        // Get patch-matching correspondence map.
        if (!std::strcmp("matchpatch",command)) {
          gmic_substitute_args(true);
          double patch_width, patch_height, patch_depth = 1, nb_iterations = 5, nb_randoms = 5;
          float patch_penalization = 0;
          unsigned int is_score = 0;
          *argx = 0; ind0.assign();
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%c",
                            gmic_use_indices,&patch_width,&end)==2 && (patch_height=patch_width)) ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%c",
                           indices,&patch_width,&patch_height,&end)==3 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf%c",
                           indices,&patch_width,&patch_height,&patch_depth,&end)==4 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf,%lf%c",
                           indices,&patch_width,&patch_height,&patch_depth,&nb_iterations,&end)==5 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf,%lf,%lf%c",
                           indices,&patch_width,&patch_height,&patch_depth,&nb_iterations,&nb_randoms,&end)==6 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf,%lf,%lf,%f%c",
                           indices,&patch_width,&patch_height,&patch_depth,&nb_iterations,&nb_randoms,
                           &patch_penalization,&end)==7 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf,%lf,%lf,%f,%u%c",
                           indices,&patch_width,&patch_height,&patch_depth,&nb_iterations,&nb_randoms,
                           &patch_penalization,&is_score,&end)==8 ||
               (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%lf,%lf,%lf,%lf,%lf,%f,%u,[%255[a-zA-Z0-9_.%+-]%c%c",
                            indices,&patch_width,&patch_height,&patch_depth,&nb_iterations,&nb_randoms,
                            &patch_penalization,&is_score,gmic_use_argx,&sep,&end)==10 && sep==']')) &&
              (ind=selection2cimg(indices,images.size(),images_names,"matchpatch")).height()==1 &&
              (!*argx ||
               (ind0=selection2cimg(argx,images.size(),images_names,"matchpatch")).height()==1) &&
              patch_width>=1 && patch_height>=1 && patch_depth>=1 &&
              nb_iterations>=0 && nb_randoms>=0 && is_score<=1) {
            const CImg<T> *initialization = 0;
            patch_width = cimg::round(patch_width);
            patch_height = cimg::round(patch_height);
            patch_depth = cimg::round(patch_depth);
            nb_iterations = cimg::round(nb_iterations);
            nb_randoms = cimg::round(nb_randoms);
            if (ind0) initialization = &images[*ind0];
            print(images,0,"Estimate correspondence map between image%s and patch image [%u], "
                  "using %gx%gx%g patches, %g iteration%s, %g randomization%s and occurrence penalization %g "
                  "(%sscore returned).",
                  gmic_selection.data(),
                  *ind,
                  patch_width,patch_height,patch_depth,
                  nb_iterations,nb_iterations!=1?"s":"",
                  nb_randoms,nb_randoms!=1?"s":"",
                  patch_penalization,
                  is_score?"":"no ");
            const CImg<T> patch_image = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply(gmic_matchpatch(patch_image,
                                                              (unsigned int)patch_width,
                                                              (unsigned int)patch_height,
                                                              (unsigned int)patch_depth,
                                                              (unsigned int)nb_iterations,
                                                              (unsigned int)nb_randoms,
                                                              patch_penalization,
                                                              (bool)is_score,
                                                              initialization),false);
          } else arg_error("matchpatch");
          is_change = true;
          ++position;
          continue;
        }

        // Project matrix onto dictionary.
        if (!std::strcmp("mproj",command)) {
          gmic_substitute_args(true);
          int method = 0, max_iter = 0;
          sep = *indices = 0; value = 1e-6;
          if ((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                           gmic_use_indices,&sep,&end)==2 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c,%d%c",
                           indices,&sep,&method,&end)==3 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c,%d,%d%c",
                           indices,&sep,&method,&max_iter,&end)==4 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c,%d,%d,%lf%c",
                           indices,&sep,&method,&max_iter,&value,&end)==5) &&
              sep==']' && method>=0 && max_iter>=0 && value>=0 &&
              (ind=selection2cimg(indices,images.size(),images_names,"mproj")).height()==1) {
            const CImg<double> A = gmic_image_arg(*ind);
            if (method==0)
              print(images,0,"Project matri%s%s to dictionary [%d] with orthogonal projection.",
                    selection.size()>1?"ce":"x",gmic_selection.data(),*ind);
            else if (method<4)
              print(images,0,"Project matri%s%s to dictionary [%d] with %s, "
                    "max iterations %d and max residual %g.",
                    selection.size()>1?"ce":"x",gmic_selection.data(),*ind,
                    method==1?"matching pursuit":
                    method==2?"matching pursuit + orthogonal projection":
                    "orthogonal matching pursuit (ortho-projection every iteration)",
                    max_iter?max_iter:A.width(),value);
            else
              print(images,0,"Project matri%s%s to dictionary [%d] with orthogonal matching pursuit "
                    "(ortho-projection every %d iterations), max iterations %d and max residual %g.",
                    selection.size()>1?"ce":"x",gmic_selection.data(),*ind,
                    method - 2,max_iter?max_iter:A.width(),value);

            cimg_forY(selection,l) gmic_apply_double(project_matrix(A,method,max_iter,value));
          } else arg_error("mproj");
          is_change = true;
          ++position;
          continue;
        }

        // Matrix division.
        gmic_arithmetic_command("mdiv",
                                operator/=,
                                "Divide matrix/vector%s by %g%s",
                                gmic_selection.data(),value,ssep,Tfloat,
                                operator/=,
                                "Divide matrix/vector%s by matrix/vector image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_diveq,
                                "Divide matrix/vector%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Divide matrix/vector%s");

        // Manage mutexes.
        if (!is_get && !std::strcmp("mutex",item)) {
          gmic_substitute_args(false);
          unsigned int number, is_lock = 1;
          if ((cimg_sscanf(argument,"%u%c",
                           &number,&end)==1 ||
               cimg_sscanf(argument,"%u,%u%c",
                           &number,&is_lock,&end)==2) &&
              number<256 && is_lock<=1) {
            print(images,0,"%s mutex #%u.",
                  is_lock?"Lock":"Unlock",number);
            if (is_lock) gmic_mutex().lock(number);
            else gmic_mutex().unlock(number);
          } else arg_error("mutex");
          ++position;
          continue;
        }

        // Draw mandelbrot/julia fractal.
        if (!std::strcmp("mandelbrot",command)) {
          gmic_substitute_args(false);
          float z0r = -2, z0i = -2, z1r = 2, z1i = 2, paramr = 0, parami = 0;
          double itermax = 100;
          unsigned int is_julia = 0;
          opacity = 1;
          if ((cimg_sscanf(argument,"%f,%f,%f,%f%c",
                           &z0r,&z0i,&z1r,&z1i,&end)==4 ||
               cimg_sscanf(argument,"%f,%f,%f,%f,%lf%c",
                           &z0r,&z0i,&z1r,&z1i,&itermax,&end)==5 ||
               cimg_sscanf(argument,"%f,%f,%f,%f,%lf,%u%c",
                           &z0r,&z0i,&z1r,&z1i,&itermax,&is_julia,&end)==6 ||
               cimg_sscanf(argument,"%f,%f,%f,%f,%lf,%u,%f,%f%c",
                           &z0r,&z0i,&z1r,&z1i,&itermax,&is_julia,&paramr,
                           &parami,&end)==8 ||
               cimg_sscanf(argument,"%f,%f,%f,%f,%lf,%u,%f,%f,%f%c",
                           &z0r,&z0i,&z1r,&z1i,&itermax,&is_julia,
                           &paramr,&parami,&opacity,&end)==9) &&
              itermax>=0 && is_julia<=1) {
            itermax = cimg::round(itermax);
            print(images,0,"Draw %s fractal on image%s, from complex area (%g,%g)-(%g,%g) "
                  "with c0 = (%g,%g) and %g iterations.",
                  is_julia?"julia":"mandelbrot",
                  gmic_selection.data(),
                  z0r,z0i,
                  z1r,z1i,
                  paramr,parami,
                  itermax);
            cimg_forY(selection,l) gmic_apply(draw_mandelbrot(CImg<T>(),opacity,z0r,z0i,z1r,z1i,(unsigned int)itermax,
                                                              true,(bool)is_julia,paramr,parami),true);
          } else arg_error("mandelbrot");
          is_change = true;
          ++position;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'n...'
        //-----------------------------
      gmic_commands_n :

        // Set image name.
        if (!is_get && !std::strcmp("name",command)) {
          gmic_substitute_args(false);
          if (is_selection && !selection)
            error(true,images,0,0,
                  "Command 'name': Empty image selection is not allowed.");

          // Extract list of specified arguments.
          const char *p = argument, *np;
          while (*p) {
            np = p; while (*np && *np!=',') ++np;
            CImg<char>(p,(unsigned int)(np - p + 1),1,1,1,true).move_to(g_list_c);
            g_list_c.back().back() = 0;
            p = np + (*np?1:0);
            if (*np && !*p) CImg<char>::vector(0).move_to(g_list_c);
          }
          if (!g_list_c) CImg<char>::vector(0).move_to(g_list_c);

          // Check correctness of specified arguments.
          if (is_selection && g_list_c.size()>selection._height)
            error(true,images,0,0,
                  "Command 'name': Number of arguments (%u) cannot be higher than "
                  "the number of images in the selection (%u).",
                  g_list_c.size(),selection._height);
          else if (!is_selection) {
            if (g_list_c.size()>images._width)
              error(true,images,0,0,
                    "Command 'name': Number of arguments (%u) cannot be higher than "
                    "the number of images in the list (%u).",
                    g_list_c.size(),images._width);

            // Retrieve implicit selection.
            selection.assign(1,g_list_c.size());
            cimg_forY(selection,l) selection[l] = images._width - g_list_c.size() + l;
          }
          print(images,0,"Set name%s of image%s to '%s'.",
                selection.height()>1?"s":"",gmic_selection.data(),gmic_argument_text_printed());
          cimglist_for(g_list_c,l) strreplace_fw(g_list_c[l]);

          cimg_forY(selection,l)
            images_names[selection[l]].assign(g_list_c[l%g_list_c.width()]);
          g_list_c.assign();
          ++position;
          continue;
        }

        // Get image indices from names.
        if (!is_get && !std::strcmp("named",command)) {
          gmic_substitute_args(false);
          if (cimg_sscanf(argument,"%u%c",&pattern,&sep)==2 && pattern<=5 && sep==',') is_cond = true;
          else { pattern = 0; is_cond = false; }
          boundary = pattern%3;
          CImg<char>::string(argument + (is_cond?2:0)).get_split(CImg<char>::vector(','),0,false).move_to(g_list_c);

          // Detect possible empty names in argument list.
          bool contains_empty_name = false;
          unsigned int sl = 0;
          for (; argument[sl]; ++sl) if (argument[sl]==',' && argument[sl + 1]==',') {
              contains_empty_name = true;
              break;
            }
          if (!contains_empty_name && sl && (*argument==',' || argument[sl - 1]==',')) contains_empty_name = true;
          if (contains_empty_name) CImg<char>(1,1,1,1,0).move_to(g_list_c);

          print(images,0,"Get %s%s with name%s '%s' for image%s (case-%s).",
                boundary==0?"all image ind":boundary==1?"lowest image ind":"highest image ind",
                boundary==0 || g_list_c.size()>1?"ices":"ex",
                g_list_c.size()>1?"s":"",
                gmic_argument_text_printed() + (is_cond?2:0),gmic_selection.data(),
                pattern<3?"sensitive":"insensitive");
          int nb_found = 0, last_found = 0;
          const bool is_single_res = boundary>0 && g_list_c.size()==1;
          if (!is_single_res) ind.assign(selection.height(),1,1,1,0);

          cimglist_for(g_list_c,k) {
            g_list_c[k].unroll('x');
            if (g_list_c[k].back()) g_list_c[k].resize(g_list_c[k].width() + 1,1,1,1,0);
            strreplace_fw(g_list_c[k]);
            switch (pattern) {
            case 0 : // All indices, case sensitive
              cimg_forY(selection,l)
                if (!std::strcmp(g_list_c[k],images_names[selection[l]])) {
                  nb_found+=(ind[l]==0); ind[l] = 1; last_found = l;
                }
              break;
            case 1 : // Lowest index, case sensitive
              if (is_single_res) {
                cimg_forY(selection,l)
                  if (!std::strcmp(g_list_c[k],images_names[selection[l]])) {
                    nb_found = 1; last_found = l; break;
                  }
              } else
                cimg_forY(selection,l)
                  if (!std::strcmp(g_list_c[k],images_names[selection[l]])) {
                    nb_found+=(ind[l]==0); ind[l] = 1; last_found = l; break;
                  }
              break;
            case 2 : // Highest index, case sensitive
              if (is_single_res) {
                cimg_rofY(selection,l)
                  if (!std::strcmp(g_list_c[k],images_names[selection[l]])) {
                    nb_found = 1; last_found = l; break;
                  }
              } else
                cimg_rofY(selection,l)
                  if (!std::strcmp(g_list_c[k],images_names[selection[l]])) {
                    nb_found+=(ind[l]==0); ind[l] = 1; last_found = l; break;
                  }
              break;
            case 3 : // All indices, case insensitive
              cimg_forY(selection,l)
                if (!cimg::strcasecmp(g_list_c[k],images_names[selection[l]])) {
                  nb_found+=(ind[l]==0); ind[l] = 1; last_found = l;
                }
              break;
            case 4 : // Lowest index, case insensitive
              if (is_single_res) {
                cimg_forY(selection,l)
                  if (!cimg::strcasecmp(g_list_c[k],images_names[selection[l]])) {
                    nb_found = 1; last_found = l; break;
                  }
              } else
                cimg_forY(selection,l)
                  if (!cimg::strcasecmp(g_list_c[k],images_names[selection[l]])) {
                    nb_found+=(ind[l]==0); ind[l] = 1; last_found = l; break;
                  }
              break;
            default : // Highest index, case insensitive
              if (is_single_res) {
                cimg_rofY(selection,l)
                  if (!cimg::strcasecmp(g_list_c[k],images_names[selection[l]])) {
                    nb_found = 1; last_found = l; break;
                  }
              } else
                cimg_rofY(selection,l)
                  if (!cimg::strcasecmp(g_list_c[k],images_names[selection[l]])) {
                    nb_found+=(ind[l]==0); ind[l] = 1; last_found = l; break;
                  }
              break;
            }
          }
          if (!nb_found) CImg<char>(1,1,1,1,0).move_to(status);
          else if (nb_found==1) {
            gmic_use_title;
            cimg_snprintf(title,_title.width(),"%u",selection[last_found]);
            CImg<char>::string(title).move_to(status);
          } else {
            ind0.assign(nb_found);
            nb_found = 0; cimg_forX(ind,l) if (ind[l]) ind0[nb_found++] = selection[l];
            ind0.value_string().move_to(status);
          }
          if (!is_single_res) ind.assign();
          g_list_c.assign();
          ++position;
          continue;
        }

        // Normalize.
        if (!std::strcmp("normalize",command)) {
          gmic_substitute_args(true);
          ind0.assign(); ind1.assign();
          sep0 = sep1 = *argx = *argy = *indices = 0;
          value0 = value1 = value = 0;
          if ((cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-]%c",
                           gmic_use_argx,gmic_use_argy,&end)==2 ||
               cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],%lf%c",
                           argx,argy,&value,&end)==3) &&
              ((cimg_sscanf(argx,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep0,&end)==2 &&
                sep0==']' &&
                (ind0=selection2cimg(indices,images.size(),images_names,"normalize")).height()==1) ||
               (cimg_sscanf(argx,"%lf%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
               cimg_sscanf(argx,"%lf%c",&value0,&end)==1) &&
              ((cimg_sscanf(argy,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_formula,&sep1,&end)==2 &&
                sep1==']' &&
                (ind1=selection2cimg(formula,images.size(),images_names,"normalize")).height()==1) ||
               (cimg_sscanf(argy,"%lf%c%c",&value1,&sep1,&end)==2 && sep1=='%') ||
               cimg_sscanf(argy,"%lf%c",&value1,&end)==1)) {
            if (ind0) { value0 = images[*ind0].min(); sep0 = 0; }
            if (ind1) { value1 = images[*ind1].max(); sep1 = 0; }
            print(images,0,"Normalize image%s in range [%g%s,%g%s], with constant-case ratio %g.",
                  gmic_selection.data(),
                  value0,sep0=='%'?"%":"",
                  value1,sep1=='%'?"%":"",
                  value);
            cimg_forY(selection,l) {
              CImg<T>& img = gmic_check(images[selection[l]]);
              nvalue0 = value0; nvalue1 = value1;
              vmin = vmax = 0;
              if (sep0=='%' || sep1=='%') {
                if (img) vmax = (double)img.max_min(vmin);
                if (sep0=='%') nvalue0 = vmin + (vmax - vmin)*value0/100;
                if (sep1=='%') nvalue1 = vmin + (vmax - vmin)*value1/100;
              }
              gmic_apply(normalize((T)nvalue0,(T)nvalue1,(float)value),true);
            }
          } else if (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep0,&end)==2 &&
                     sep0==']' &&
                     (ind0=selection2cimg(indices,images.size(),images_names,"normalize")).height()==1) {
            if (images[*ind0]) value1 = (double)images[*ind0].max_min(value0);
            print(images,0,"Normalize image%s in range [%g,%g].",
                  gmic_selection.data(),
                  value0,
                  value1);
            cimg_forY(selection,l) gmic_apply(normalize((T)value0,(T)value1),true);
          } else arg_error("normalize");
          is_change = true;
          ++position;
          continue;
        }

        // Test difference.
        gmic_arithmetic_command("neq",
                                operator_neq,
                                "Compute boolean inequality between image%s and %g%s",
                                gmic_selection.data(),value,ssep,T,
                                operator_neq,
                                "Compute boolean inequality between image%s and image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_neq,
                                "Compute boolean inequality between image%s and expression %s'",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute boolean inequality between image%s");

        // Manage network permission and timeout.
        if (!is_get && !std::strcmp("network",item)) {
          gmic_substitute_args(false);
          if (cimg_sscanf(argument,"%d%c",
                          &err,&end)==1 &&
              err>=-1) {
            if (err==-1) print(images,0,"Disable load-from-network.");
            else if (!err) print(images,0,"Enable load-from-network, with no timeout.");
            else print(images,0,"Enable load-from-network, with %ds timeout.",err);
            cimg::network_mode(err!=-1,true);
            if (err!=-1) network_timeout = err;
          } else arg_error("network");
          ++position;
          continue;
        }

        // Discard custom command arguments.
        if (!is_get && !std::strcmp("noarg",item)) {
          print(images,0,"Discard command arguments.");
          if (is_noarg) *is_noarg = true;
          continue;
        }

        // Add noise.
        if (!std::strcmp("noise",command)) {
          gmic_substitute_args(false);
          int noise_type = 0;
          float sigma = 0;
          sep = 0;
          if ((cimg_sscanf(argument,"%f%c",
                           &sigma,&end)==1 ||
               (cimg_sscanf(argument,"%f%c%c",
                            &sigma,&sep,&end)==2 && sep=='%') ||
               cimg_sscanf(argument,"%f,%d%c",
                           &sigma,&noise_type,&end)==2 ||
               (cimg_sscanf(argument,"%f%c,%d%c",
                            &sigma,&sep,&noise_type,&end)==3 && sep=='%')) &&
              sigma>=0 && noise_type>=0 && noise_type<=4) {
            const char *s_type = noise_type==0?"gaussian":
              noise_type==1?"uniform":
              noise_type==2?"salt&pepper":
              noise_type==3?"poisson":"rice";
            if (sep=='%' && noise_type!=2) sigma = -sigma;
            print(images,0,"Add %s noise to image%s, with standard deviation %g%s.",
                  s_type,
                  gmic_selection.data(),
                  cimg::abs(sigma),sep=='%'?"%":"");
            cimg_forY(selection,l) gmic_apply(noise(sigma,noise_type),true);
          } else arg_error("noise");
          is_change = true;
          ++position;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'o...'
        //-----------------------------
      gmic_commands_o :

        // Exception handling in local environments.
        if (!is_get && !std::strcmp("onfail",item)) {
          const CImg<char> &s = callstack.back();
          if (s[0]!='*' || s[1]!='l')
            error(true,images,0,0,
                  "Command 'onfail': Not associated to a 'local' command within the same scope.");
          for (int nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
            it = commands_line[position];
            if (*it==1)
              is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
            else {
              _is_get = *it=='+';
              it+=(_is_get || *it=='-');
              gmic_if_flr ++nb_levels; gmic_elif_flr { if (!--nb_levels) --position; }
            }
          }
          continue;
        }

        // Draw 3D object.
        if (!std::strcmp("object3d",command)) {
          gmic_substitute_args(true);
          unsigned int is_zbuffer = 1, _double3d = ~0U, _render3d = ~0U;
          float x = 0, y = 0, z = 0,
            _focale3d = cimg::type<float>::nan(),
            _specl3d = cimg::type<float>::nan(),
            _specs3d = cimg::type<float>::nan(),
            _light3d_x = light3d_x,
            _light3d_y = light3d_y,
            _light3d_z = light3d_z;
          sep = sepx = sepy = *argx = *argy = 0;
          opacity = 1;
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                            gmic_use_indices,&sep,&end)==2 && sep==']') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-]%c",
                           indices,gmic_use_argx,&end)==2 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           indices,argx,gmic_use_argy,&end)==3 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%f%c",
                           indices,argx,argy,&z,&end)==4 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%f,%f%c",
                           indices,argx,argy,&z,&opacity,&end)==5 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%f,%f,%u%c",
                           indices,argx,argy,&z,&opacity,&_render3d,&end)==6 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%f,%f,%u,%u%c",
                           indices,argx,argy,&z,&opacity,&_render3d,&_double3d,&end)==7 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%f,%f,%u,%u,%u%c",
                           indices,argx,argy,&z,&opacity,&_render3d,&_double3d,&is_zbuffer,
                           &end)==8 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%f,%f,%u,%u,%u,%f%c",
                           indices,argx,argy,&z,&opacity,&_render3d,&_double3d,&is_zbuffer,
                           &_focale3d,&end)==9 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%f,%f,%u,%u,%u,%f,%f,%f,%f%c",
                           indices,argx,argy,&z,&opacity,&_render3d,&_double3d,&is_zbuffer,
                           &_focale3d,&_light3d_x,&_light3d_y,&_light3d_z,&end)==12 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%f,%f,%u,%u,%u,%f,%f,%f,%f,%f%c",
                           indices,argx,argy,&z,&opacity,&_render3d,&_double3d,&is_zbuffer,
                           &_focale3d,&_light3d_x,&_light3d_y,&_light3d_z,
                           &_specl3d,&end)==13 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%f,%f,%u,%u,%u,%f,%f,%f,%f,%f,%f%c",
                           indices,argx,argy,&z,&opacity,&_render3d,&_double3d,&is_zbuffer,
                           &_focale3d,&_light3d_x,&_light3d_y,&_light3d_z,
                           &_specl3d,&_specs3d,&end)==14) &&
              (ind=selection2cimg(indices,images.size(),images_names,"object3d")).height()==1 &&
              (!*argx ||
               cimg_sscanf(argx,"%f%c",&x,&end)==1 ||
               (cimg_sscanf(argx,"%f%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy ||
               cimg_sscanf(argy,"%f%c",&y,&end)==1 ||
               (cimg_sscanf(argy,"%f%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
              (_render3d==~0U || _render3d<=5) && is_zbuffer<=1 &&
              (_double3d==~0U || _double3d<=1)) {

            // Get default rendering options.
            if (_render3d==~0U) { // Rendering mode
              _render3d = 4;
              const CImg<char> s_mode3d = get_variable("_mode3d",variables_sizes,0,0);
              if (s_mode3d && *s_mode3d) {
                if (*s_mode3d>='0' && *s_mode3d<='5' && !s_mode3d[1]) _render3d = (unsigned int)(*s_mode3d - '0');
                else if (*s_mode3d=='-' && s_mode3d[1]=='1' && !s_mode3d[2]) _render3d = 0;
              }
            }
            if (_double3d==~0U) { // Simple/double face property
              _double3d = 1;
              const CImg<char> s_double3d = get_variable("_double3d",variables_sizes,0,0);
              if (s_double3d && *s_double3d>='0' && *s_double3d<='1' && !s_double3d[1])
                _double3d = (unsigned int)(*s_double3d - '0');
            }
            if (cimg::type<float>::is_nan(_focale3d)) { // Focale
              const CImg<char> s_focale3d = get_variable("_focale3d",variables_sizes,0,0);
              if (!s_focale3d || std::sscanf(s_focale3d,"%f%c",&_focale3d,&end)!=1)
                _focale3d = 700;
            }
            if (cimg::type<float>::is_nan(_specl3d)) { // Specular lightness
              const CImg<char> s_specl3d = get_variable("_specl3d",variables_sizes,0,0);
              if (!s_specl3d || std::sscanf(s_specl3d,"%f%c",&_specl3d,&end)!=1 || _specl3d<0)
                _specl3d = 0.15f;
            }
            if (cimg::type<float>::is_nan(_specs3d)) { // Specular shininess
              const CImg<char> s_specs3d = get_variable("_specs3d",variables_sizes,0,0);
              if (!s_specs3d || std::sscanf(s_specs3d,"%f%c",&_specs3d,&end)!=1 || _specs3d<0)
                _specs3d = 0.15f;
            }

            print(images,0,"Draw 3D object [%u] at (%g%s,%g%s,%g) on image%s, with opacity %g, "
                  "%s rendering, %s-sided mode, %sz-buffer, focale %g, 3D light at (%g,%g,%g) "
                  "and specular properties (%g,%g)",
                  *ind,
                  x,sepx=='%'?"%":"",
                  y,sepy=='%'?"%":"",
                  z,
                  gmic_selection.data(),
                  opacity,
                  _render3d==0?"dot":_render3d==1?"wireframe":_render3d==2?"flat":
                  _render3d==3?"flat-shaded":_render3d==4?"gouraud-shaded":"phong-shaded",
                  _double3d?"double":"simple",
                  is_zbuffer?"":"no ",
                  _focale3d,_light3d_x,_light3d_y,_light3d_z,
                  _specl3d,_specs3d);

            const CImg<T> img0 = gmic_image_arg(*ind);
            CImgList<float> opacities;
            vertices.assign(img0,false);

            try {
              if (_render3d>=3) {
                vertices.CImg3dtoobject3d(primitives,g_list_uc,opacities,false);
                if (light3d) g_list_uc.insert(light3d,~0U,true);
              } else vertices.CImg3dtoobject3d(primitives,g_list_f,opacities,false);
            } catch (CImgException&) {
              if (!vertices.is_CImg3d(true,&(*gmic_use_message=0)))
                error(true,images,0,0,
                      "Command 'object3d': Invalid 3D object [%u], specified "
                      "in argument '%s' (%s).",
                      *ind,gmic_argument_text(),message);
              else throw;
            }

            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const float
                nx = sepx=='%'?x*(img.width() - 1)/100:x,
                ny = sepy=='%'?y*(img.height() - 1)/100:y;
              CImg<float> zbuffer(is_zbuffer?img.width():0,is_zbuffer?img.height():0,1,1,0);
              if (g_list_f) {
                gmic_apply(draw_object3d(nx,ny,z,vertices,primitives,g_list_f,opacities,
                                         _render3d,_double3d,_focale3d,
                                         _light3d_x,_light3d_y,_light3d_z,
                                         _specl3d,_specs3d,
                                         opacity,zbuffer),true);

              } else {
                gmic_apply(draw_object3d(nx,ny,z,vertices,primitives,g_list_uc,opacities,
                                         _render3d,_double3d,_focale3d,
                                         _light3d_x,_light3d_y,_light3d_z,
                                         _specl3d,_specs3d,
                                         opacity,zbuffer),true);
              }
            }
          } else arg_error("object3d");
          g_list_f.assign();
          g_list_uc.assign();
          vertices.assign();
          primitives.assign();
          is_change = true;
          ++position;
          continue;
        }

        // Bitwise or.
        gmic_arithmetic_command("or",
                                operator|=,
                                "Compute bitwise OR of image%s by %g%s",
                                gmic_selection.data(),value,ssep,Tlong,
                                operator|=,
                                "Compute bitwise OR of image%s by image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_oreq,
                                "Compute bitwise OR of image%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute sequential bitwise OR of image%s");

        // Output.
        if (!is_get && !std::strcmp("output",command)) {
          gmic_substitute_args(false);

          // Set good alias for shared variables.
          gmic_use_argc;
          gmic_use_title;
          gmic_use_color;
          CImg<char> &_filename = _color, &filename_tmp = _title, &options = _argc;
          char cext[12];
          *cext = *_filename = *filename_tmp = *options = 0;
          CImgList<unsigned int> empty_indices;
          CImg<char> eselec;

          if (cimg_sscanf(argument,"%11[a-zA-Z0-9]:%4095[^,],%255s", // Detect forced file format
                          cext,_filename.data(),options.data())<2 ||
              !cext[1]) { // Length of preprend 'ext' must be >=2 (avoid case 'C:\\...' on Windows)
            *cext = *_filename = *options = 0;
            if (cimg_sscanf(argument,"%4095[^,],%255s",_filename.data(),options.data())!=2) {
              std::strncpy(_filename,argument,_filename.width() - 1);
              _filename[_filename.width() - 1] = 0;
            }
          }
          strreplace_fw(_filename);
          strreplace_fw(options);
          const bool is_stdout = *_filename=='-' && (!_filename[1] || _filename[1]=='.');

          if (*cext) { // Force output to be written as a '.ext' file : generate random filename
            if (is_stdout) {
              // Simplify filename 'ext:-.foo' as '-.ext'.
              cimg_snprintf(_filename,_filename.width(),"-.%s",cext);
              *cext = 0;
            } else {
              std::FILE *file = 0;
              do {
                cimg_snprintf(filename_tmp,filename_tmp.width(),"%s%c%s.%s",
                              cimg::temporary_path(),cimg_file_separator,
                              cimg::filenamerand(),cext);
                if ((file=cimg::std_fopen(filename_tmp,"rb"))!=0) cimg::fclose(file);
              } while (file);
            }
          }
          const char
            *const filename = *cext?filename_tmp:_filename,
            *const ext = cimg::split_filename(filename);
          CImg<char> uext = CImg<char>::string(ext);
          cimg::lowercase(uext);

          if (!std::strcmp(uext,"off")) {

            // OFF file (geomview).
            *gmic_use_formula = 0;
            std::strncpy(formula,filename,_formula.width() - 1);
            formula[_formula.width() - 1] = 0;

            if (*options)
              error(true,images,0,0,
                    "Command 'output': File '%s', format does not take any output options (options '%s' specified).",
                    formula,options.data());

            print(images,0,"Output 3D object%s as %s file '%s'.",
                  gmic_selection.data(),uext.data(),formula);

            cimg_forY(selection,l) {
              uind = selection[l];
              const CImg<T>& img = gmic_check(images[uind]);
              if (selection.height()!=1) cimg::number_filename(filename,l,6,formula);
              CImgList<float> opacities;
              vertices.assign(img,false);
              try {
                vertices.CImg3dtoobject3d(primitives,g_list_f,opacities,false).
                  save_off(primitives,g_list_f,formula);
              } catch (CImgException&) {
                if (!vertices.is_CImg3d(true,&(*gmic_use_message=0)))
                  error(true,images,0,0,
                        "Command 'output': 3D object file '%s', invalid 3D object [%u] "
                        "in image%s (%s).",
                        formula,uind,gmic_selection.data(),message);
                else throw;
              }
            }
            vertices.assign();
            primitives.assign();
            g_list_f.assign();
          } else if (!std::strcmp(uext,"cpp") || !std::strcmp(uext,"c") ||
                     !std::strcmp(uext,"hpp") || !std::strcmp(uext,"h") ||
                     !std::strcmp(uext,"pan") || !std::strcmp(uext,"pnk") ||
                     !std::strcmp(uext,"pgm") || !std::strcmp(uext,"ppm") ||
                     !std::strcmp(uext,"ppm") || !std::strcmp(uext,"hdr") ||
                     !std::strcmp(uext,"nii")) {

            // .cpp, .c, .hpp, .h, .pan, .pnk, .pgm, .ppm, .pnm, .hdr or .nii file.
            const char *
              stype = (cimg_sscanf(options,"%255[a-z123468]%c",&(*gmic_use_argx=0),&(end=0))==1 ||
                       (cimg_sscanf(options,"%255[a-z123468]%c",&(*argx=0),&end)==2 && end==','))?
              argx:"auto";
            g_list.assign(selection.height());
            cimg_forY(selection,l) if (!gmic_check(images[selection(l)]))
              CImg<unsigned int>::vector(selection(l)).move_to(empty_indices);
            if (empty_indices && is_verbose) {
              selection2string(empty_indices>'y',images_names,1,eselec);
              warn(images,0,false,
                   "Command 'output': Image%s %s empty.",
                   eselec.data(),empty_indices.size()>1?"are":"is");
            }
            cimg_forY(selection,l)
              g_list[l].assign(images[selection[l]],true);
            if (g_list.size()==1)
              print(images,0,
                    "Output image%s as %s file '%s', with pixel type '%s' (1 image %dx%dx%dx%d).",
                    gmic_selection.data(),
                    uext.data(),_filename.data(),
                    stype,
                    g_list[0].width(),g_list[0].height(),
                    g_list[0].depth(),g_list[0].spectrum());
            else print(images,0,"Output image%s as %s file '%s', with pixel type '%s'.",
                       gmic_selection.data(),
                       uext.data(),_filename.data(),
                       stype);
            if (!g_list)
              error(true,images,0,0,
                    "Command 'output': File '%s', instance list (%u,%p) is empty.",
                    _filename.data(),g_list.size(),g_list.data());

#define gmic_save_multitype(svalue_type,value_type) \
              if (!std::strcmp(stype,svalue_type)) { \
                if (g_list.size()==1) \
                  CImg<value_type>::copy_rounded(g_list[0]).save(filename); \
                else { \
                  cimglist_for(g_list,l) { \
                    cimg::number_filename(filename,l,6,gmic_use_formula); \
                    CImg<value_type>::copy_rounded(g_list[l]).save(formula); \
                  } \
                } \
              }
            if (!std::strcmp(stype,"auto")) stype = CImg<T>::storage_type(g_list,false);
            gmic_save_multitype("uint8",cimg_uint8)
            else gmic_save_multitype("int8",cimg_int8)
              else gmic_save_multitype("uint16",cimg_uint16)
                else gmic_save_multitype("int16",cimg_int16)
                  else gmic_save_multitype("uint32",cimg_uint32)
                    else gmic_save_multitype("int32",cimg_int32)
                      else gmic_save_multitype("uint64",cimg_uint64)
                        else gmic_save_multitype("int64",cimg_int64)
                          else gmic_save_multitype("float32",cimg_float32)
                            else gmic_save_multitype("float64",cimg_float64)
                              else error(true,images,0,0,
                                         "Command 'output': File '%s', invalid "
                                         "specified pixel type '%s'.",
                                         _filename.data(),stype);
          } else if (!std::strcmp(uext,"tiff") || !std::strcmp(uext,"tif")) {

            // TIFF file.
            const char *
              stype = (cimg_sscanf(options,"%255[a-z123468]%c",&(*gmic_use_argx=0),&(end=0))==1 ||
                       (cimg_sscanf(options,"%255[a-z123468]%c",&(*argx=0),&end)==2 && end==','))?
              argx:"auto";
            const unsigned int l_stype = (unsigned int)std::strlen(stype);
            const char *const _options = options.data() + (stype!=argx?0:l_stype + (end==','?1:0));
            float _is_multipage = 0;
            *argy = 0; opacity = 1;
            if (cimg_sscanf(_options,"%255[a-z],%f,%f",gmic_use_argy,&_is_multipage,&opacity)<1)
              cimg_sscanf(_options,"%f,%f",&_is_multipage,&opacity);
            const unsigned int compression_type =
              !cimg::strcasecmp(argy,"jpeg") ||
              !cimg::strcasecmp(argy,"jpg")?2:
              !cimg::strcasecmp(argy,"lzw")?1U:0U;
            const bool is_multipage = (bool)cimg::round(_is_multipage);
            const bool use_bigtiff = (bool)cimg::round(opacity);

            g_list.assign(selection.height());
            cimg_forY(selection,l) if (!gmic_check(images[selection(l)]))
              CImg<unsigned int>::vector(selection(l)).move_to(empty_indices);
            if (empty_indices && is_verbose) {
              selection2string(empty_indices>'y',images_names,1,eselec);
              warn(images,0,false,
                   "Command 'output': Image%s %s empty.",
                   eselec.data(),empty_indices.size()>1?"are":"is");
            }
            cimg_forY(selection,l)
              g_list[l].assign(images[selection[l]],true);
            if (g_list.size()==1)
              print(images,0,"Output image%s as %s file '%s', with pixel type '%s', %s compression "
                    "and %sbigtiff support (1 image %dx%dx%dx%d).",
                    gmic_selection.data(),
                    uext.data(),_filename.data(),stype,
                    compression_type==2?"JPEG":compression_type==1?"LZW":"no",
                    use_bigtiff?"":"no ",
                    g_list[0].width(),g_list[0].height(),
                    g_list[0].depth(),g_list[0].spectrum());
            else print(images,0,"Output image%s as %s file '%s', with pixel type '%s', "
                       "%s compression, %s-page mode and %s bigtiff support.",
                       gmic_selection.data(),
                       uext.data(),_filename.data(),stype,
                       compression_type==2?"JPEG":compression_type==1?"LZW":"no",
                       is_multipage?"multi":"single",
                       use_bigtiff?"":"no ");
            if (!g_list)
              error(true,images,0,0,
                    "Command 'output': File '%s', instance list (%u,%p) is empty.",
                    _filename.data(),g_list.size(),g_list.data());

#define gmic_save_tiff(svalue_type,value_type) \
              if (!std::strcmp(stype,svalue_type)) { \
                if (g_list.size()==1 || is_multipage) \
                  CImgList<value_type>::copy_rounded(g_list). \
                    save_tiff(filename,compression_type,0,0,use_bigtiff); \
                else { \
                  cimglist_for(g_list,l) { \
                    cimg::number_filename(filename,l,6,gmic_use_formula); \
                    CImg<value_type>::copy_rounded(g_list[l]). \
                      save_tiff(formula,compression_type,0,0,use_bigtiff); \
                  } \
                } \
              }
            if (!std::strcmp(stype,"auto")) stype = CImg<T>::storage_type(g_list,false);
            gmic_save_tiff("uint8",cimg_uint8)
            else gmic_save_tiff("int8",cimg_int8)
              else gmic_save_tiff("uint16",cimg_uint16)
                else gmic_save_tiff("int16",cimg_int16)
                  else gmic_save_tiff("uint32",cimg_uint32)
                    else gmic_save_tiff("int32",cimg_int32)
                      else gmic_save_tiff("uint64",cimg_uint64)
                        else gmic_save_tiff("int64",cimg_int64)
                          else gmic_save_tiff("float32",cimg_float32)
                            else gmic_save_tiff("float64",cimg_float64)
                              else error(true,images,0,0,
                                         "Command 'output': File '%s', invalid "
                                         "specified pixel type '%s'.",
                                         _filename.data(),stype);

          } else if (!std::strcmp(uext,"gif")) {

            // GIF file.
            float fps = 0, _nb_loops = 0;
            g_list.assign(selection.height());
            cimg_forY(selection,l) if (!gmic_check(images[selection(l)]))
              CImg<unsigned int>::vector(selection(l)).move_to(empty_indices);
            if (empty_indices && is_verbose) {
              selection2string(empty_indices>'y',images_names,1,eselec);
              warn(images,0,false,
                   "Command 'output': Image%s %s empty.",
                   eselec.data(),empty_indices.size()>1?"are":"is");
            }
            cimg_forY(selection,l)
              g_list[l].assign(images[selection[l]],true);
            if (g_list.size()>1 && cimg_sscanf(options,"%f,%f",&fps,&_nb_loops)>=1 && fps>0) {
              // Save animated .gif file.
              const unsigned int nb_loops = (unsigned int)cimg::round(_nb_loops);
              if (nb_loops)
                print(images,0,
                      "Output image%s as animated %s file '%s', with %g fps and %u loops.",
                      gmic_selection.data(),uext.data(),_filename.data(),fps,nb_loops);
              else
                print(images,0,
                      "Output image%s as animated %s file '%s', with %g fps.",
                      gmic_selection.data(),uext.data(),_filename.data(),fps);
              g_list.save_gif_external(filename,fps,nb_loops);
            } else {
              if (g_list.size()==1)
                print(images,0,"Output image%s as %s file '%s' (1 image %dx%dx%dx%d).",
                      gmic_selection.data(),
                      uext.data(),_filename.data(),
                      g_list[0].width(),g_list[0].height(),
                      g_list[0].depth(),g_list[0].spectrum());
              else print(images,0,"Output image%s as %s file '%s'.",
                         gmic_selection.data(),uext.data(),_filename.data());
              g_list.save(filename); // Save distinct .gif files
            }
          } else if (!std::strcmp(uext,"jpeg") || !std::strcmp(uext,"jpg")) {

            // JPEG file.
            float quality = 100;
            if (cimg_sscanf(options,"%f%c",&quality,&end)!=1) quality = 100;
            if (quality<0) quality = 0; else if (quality>100) quality = 100;
            g_list.assign(selection.height());
            cimg_forY(selection,l) if (!gmic_check(images[selection(l)]))
              CImg<unsigned int>::vector(selection(l)).move_to(empty_indices);
            if (empty_indices && is_verbose) {
              selection2string(empty_indices>'y',images_names,1,eselec);
              warn(images,0,false,
                   "Command 'output': Image%s %s empty.",
                   eselec.data(),empty_indices.size()>1?"are":"is");
            }
            cimg_forY(selection,l)
              g_list[l].assign(images[selection[l]],true);
            if (g_list.size()==1)
              print(images,0,
                    "Output image%s as %s file '%s', with quality %g%% (1 image %dx%dx%dx%d).",
                    gmic_selection.data(),
                    uext.data(),_filename.data(),
                    quality,
                    g_list[0].width(),g_list[0].height(),
                    g_list[0].depth(),g_list[0].spectrum());
            else print(images,0,"Output image%s as %s file '%s', with quality %g%%.",
                       gmic_selection.data(),
                       uext.data(),_filename.data(),
                       quality);
            if (!g_list)
              error(true,images,0,0,
                    "Command 'output': File '%s', instance list (%u,%p) is empty.",
                    _filename.data(),g_list.size(),g_list.data());
            if (g_list.size()==1)
              g_list[0].save_jpeg(filename,(unsigned int)cimg::round(quality));
            else {
              cimglist_for(g_list,l) {
                cimg::number_filename(filename,l,6,gmic_use_formula);
                g_list[l].save_jpeg(formula,(unsigned int)cimg::round(quality));
              }
            }
          } else if (!std::strcmp(uext,"mnc") && *options) {

            // MNC file.
            g_list.assign(selection.height());
            cimg_forY(selection,l) if (!gmic_check(images[selection(l)]))
              CImg<unsigned int>::vector(selection(l)).move_to(empty_indices);
            if (empty_indices && is_verbose) {
              selection2string(empty_indices>'y',images_names,1,eselec);
              warn(images,0,false,
                   "Command 'output': Image%s %s empty.",
                   eselec.data(),empty_indices.size()>1?"are":"is");
            }
            cimg_forY(selection,l)
              g_list[l].assign(images[selection[l]],true);
            if (g_list.size()==1)
              print(images,0,
                    "Output image%s as %s file '%s', with header get from file '%s' "
                    "(1 image %dx%dx%dx%d).",
                    gmic_selection.data(),
                    uext.data(),_filename.data(),
                    options.data(),
                    g_list[0].width(),g_list[0].height(),
                    g_list[0].depth(),g_list[0].spectrum());
            else
              print(images,0,
                    "Output image%s as %s file '%s', with header get from file '%s'.",
                    gmic_selection.data(),
                    uext.data(),_filename.data(),
                    options.data());
            if (g_list.size()==1)
              g_list[0].save_minc2(filename,options);
            else {
              cimglist_for(g_list,l) {
                cimg::number_filename(filename,l,6,gmic_use_formula);
                g_list[l].save_minc2(formula,options);
              }
            }
          } else if (!std::strcmp(uext,"raw")) {

            // RAW data file.
            const char *stype = cimg_sscanf(options,"%255[a-z123468]%c",gmic_use_argx,&end)==1?argx:"auto";
            g_list.assign(selection.height());
            cimg_forY(selection,l) if (!gmic_check(images[selection(l)]))
              CImg<unsigned int>::vector(selection(l)).move_to(empty_indices);
            if (empty_indices && is_verbose) {
              selection2string(empty_indices>'y',images_names,1,eselec);
              warn(images,0,false,
                   "Command 'output': Image%s %s empty.",
                   eselec.data(),empty_indices.size()>1?"are":"is");
            }
            cimg_forY(selection,l)
              g_list[l].assign(images[selection[l]],true);
            if (g_list.size()==1)
              print(images,0,
                    "Output image%s as %s file '%s', with pixel type '%s' (1 image %dx%dx%dx%d).",
                    gmic_selection.data(),
                    uext.data(),_filename.data(),
                    stype,
                    g_list[0].width(),g_list[0].height(),
                    g_list[0].depth(),g_list[0].spectrum());
            else print(images,0,"Output image%s as %s file '%s', with pixel type '%s'.",
                       gmic_selection.data(),
                       uext.data(),_filename.data(),
                       stype);
            if (!g_list)
              error(true,images,0,0,
                    "Command 'output': File '%s', instance list (%u,%p) is empty.",
                    _filename.data(),g_list.size(),g_list.data());

#define gmic_save_raw(svalue_type,value_type) \
              if (!std::strcmp(stype,svalue_type)) { \
                if (g_list.size()==1) \
                  CImg<value_type>::copy_rounded(g_list[0]).save_raw(filename); \
                else { \
                  cimglist_for(g_list,l) { \
                    cimg::number_filename(filename,l,6,gmic_use_formula); \
                    CImg<value_type>::copy_rounded(g_list[l]).save_raw(formula); \
                  } \
                } \
              }
            if (!std::strcmp(stype,"auto")) stype = CImg<T>::storage_type(g_list,true);
            gmic_save_raw("bool",bool)
            else gmic_save_raw("uint8",cimg_uint8)
              else gmic_save_raw("int8",cimg_int8)
                else gmic_save_raw("uint16",cimg_uint16)
                  else gmic_save_raw("int16",cimg_int16)
                    else gmic_save_raw("uint32",cimg_uint32)
                      else gmic_save_raw("int32",cimg_int32)
                        else gmic_save_raw("uint64",cimg_uint64)
                          else gmic_save_raw("int64",cimg_int64)
                            else gmic_save_raw("float32",cimg_float32)
                              else gmic_save_raw("float64",cimg_float64)
                                else error(true,images,0,0,
                                           "Command 'output': File '%s', invalid "
                                           "specified pixel type '%s'.",
                                           _filename.data(),stype);
          } else if (!std::strcmp(uext,"yuv")) {

            // YUV sequence.
            if (cimg_sscanf(options,"%f",&opacity)!=1) opacity = 444;
            else opacity = cimg::round(opacity);
            const unsigned int ich = (unsigned int)opacity;
            if (ich!=420 && ich!=422 && ich!=444)
              error(true,images,0,0,
                    "Command 'output': YUV file '%s', specified chroma subsampling %g is invalid.",
                    _filename.data(),opacity);
            g_list.assign(selection.height());
            cimg_forY(selection,l)
              g_list[l].assign(images[selection[l]],true);
            print(images,0,"Output image%s as YUV-%u:%u:%u file '%s'.",
                  gmic_selection.data(),
                  ich/100,(ich/10)%10,ich%10,
                  _filename.data());
            g_list.save_yuv(filename,ich,true);

          } else if (!std::strcmp(uext,"cimg") || !std::strcmp(uext,"cimgz")) {

            // CImg[z] file.
            const char *stype = cimg_sscanf(options,"%255[a-z123468]%c",gmic_use_argx,&end)==1?argx:"auto";
            g_list.assign(selection.height());
            cimg_forY(selection,l)
              g_list[l].assign(images[selection[l]],true);
            print(images,0,"Output image%s as %s file '%s', with pixel type '%s'.",
                  gmic_selection.data(),
                  uext.data(),_filename.data(),
                  stype);

#define gmic_save_cimg(svalue_type,value_type) \
              if (!std::strcmp(stype,svalue_type)) \
                CImgList<value_type>::copy_rounded(g_list).save(filename);

            if (!std::strcmp(stype,"auto")) stype = CImg<T>::storage_type(g_list,false);
            gmic_save_cimg("uint8",cimg_uint8)
            else gmic_save_cimg("int8",cimg_int8)
              else gmic_save_cimg("uint16",cimg_uint16)
                else gmic_save_cimg("int16",cimg_int16)
                  else gmic_save_cimg("uint32",cimg_uint32)
                    else gmic_save_cimg("int32",cimg_int32)
                      else gmic_save_cimg("uint64",cimg_uint64)
                        else gmic_save_cimg("int64",cimg_int64)
                          else gmic_save_cimg("float32",cimg_float32)
                            else gmic_save_cimg("float64",cimg_float64)
                              else error(true,images,0,0,
                                         "Command 'output': File '%s', invalid "
                                         "specified pixel type '%s'.",
                                         _filename.data(),stype);
          } else if (!std::strcmp(uext,"gmz") || !*ext) {

            // GMZ file.
            const char *stype = cimg_sscanf(options,"%255[a-z123468]%c",gmic_use_argx,&end)==1?argx:"auto";
            g_list.assign(selection.height());
            g_list_c.assign(selection.height());
            cimg_forY(selection,l) {
              g_list[l].assign(images[selection[l]],true);
              g_list_c[l].assign(images_names[selection[l]],true);
            }
            print(images,0,"Output image%s as %s file '%s', with pixel type '%s'.",
                  gmic_selection.data(),
                  uext.data(),_filename.data(),
                  stype);

#define gmic_save_gmz(svalue_type,value_type) \
              if (!std::strcmp(stype,svalue_type)) \
                CImg<value_type>::save_gmz(filename,CImgList<value_type>::copy_rounded(g_list),g_list_c);
            if (!std::strcmp(stype,"auto")) stype = CImg<T>::storage_type(g_list,false);
            gmic_save_gmz("uint8",cimg_uint8)
            else gmic_save_gmz("int8",cimg_int8)
              else gmic_save_gmz("uint16",cimg_uint16)
                else gmic_save_gmz("int16",cimg_int16)
                  else gmic_save_gmz("uint32",cimg_uint32)
                    else gmic_save_gmz("int32",cimg_int32)
                      else gmic_save_gmz("uint64",cimg_uint64)
                        else gmic_save_gmz("int64",cimg_int64)
                          else gmic_save_gmz("float32",cimg_float32)
                            else gmic_save_gmz("float64",cimg_float64)
                              else error(true,images,0,0,
                                         "Command 'output': File '%s', invalid "
                                         "specified pixel type '%s'.",
                                         _filename.data(),stype);
          } else if (!std::strcmp(uext,"avi") ||
                     !std::strcmp(uext,"mov") ||
                     !std::strcmp(uext,"asf") ||
                     !std::strcmp(uext,"divx") ||
                     !std::strcmp(uext,"flv") ||
                     !std::strcmp(uext,"mpg") ||
                     !std::strcmp(uext,"m1v") ||
                     !std::strcmp(uext,"m2v") ||
                     !std::strcmp(uext,"m4v") ||
                     !std::strcmp(uext,"mjp") ||
                     !std::strcmp(uext,"mp4") ||
                     !std::strcmp(uext,"mkv") ||
                     !std::strcmp(uext,"mpe") ||
                     !std::strcmp(uext,"movie") ||
                     !std::strcmp(uext,"ogm") ||
                     !std::strcmp(uext,"ogg") ||
                     !std::strcmp(uext,"ogv") ||
                     !std::strcmp(uext,"qt") ||
                     !std::strcmp(uext,"rm") ||
                     !std::strcmp(uext,"vob") ||
                     !std::strcmp(uext,"webm") ||
                     !std::strcmp(uext,"wmv") ||
                     !std::strcmp(uext,"xvid") ||
                     !std::strcmp(uext,"mpeg")) {

            // Generic video file.
            float fps = 0, keep_open = 0;
            name.assign(8); *name = 0; // Codec
            cimg_sscanf(options,"%f,%7[a-zA-Z0-9],%f",&fps,name.data(),&keep_open);
            fps = cimg::round(fps);
            if (!fps) fps = 25;
            if (*name=='0' && !name[1]) *name = 0;
            g_list.assign(selection.height());
            cimg_forY(selection,l) if (!gmic_check(images[selection(l)]))
              CImg<unsigned int>::vector(selection(l)).move_to(empty_indices);
            if (empty_indices && is_verbose) {
              selection2string(empty_indices>'y',images_names,1,eselec);
              warn(images,0,false,
                   "Command 'output': Image%s %s empty.",
                   eselec.data(),empty_indices.size()>1?"are":"is");
            }
            cimg_forY(selection,l)
              g_list[l].assign(images[selection[l]],true);
            print(images,0,"Output image%s as %s file '%s', with %g fps and %s codec.",
                  gmic_selection.data(),
                  uext.data(),_filename.data(),
                  fps,*name?name.data():"(default)");

#ifndef cimg_use_opencv
            if (keep_open)
              warn(images,0,false,
                   "Command 'output': Cannot output streamed video, as this requires features from the "
                   "OpenCV library (not enabled at compilation time).");
#endif
            g_list.save_video(filename,(unsigned int)fps,name,(bool)keep_open);
            if (!cimg::fsize(filename)) throw CImgException("Output file '%s' is empty. Something went wrong!",
                                                            _filename.data());
          } else { // Any other extension

            // Check if a custom command handling requested file format exists.
            gmic_use_formula;
            cimg_snprintf(formula,_formula.width(),"output_%s",uext.data());
            hash = hashcode(formula,false);
            if (search_sorted(formula,commands_names[hash],commands_names[hash].size(),pattern)) { // Command found
              cimg_snprintf(formula,_formula.width(),"output_%s[%s] \"%s\"%s%s",
                            uext.data(),*s_selection?s_selection:"^",filename,
                            *options?",":"",options.data());
              const CImgList<char> ncommands_line = commands_line_to_CImgList(formula);
              unsigned int nposition = 0;
              CImg<char>::string("").move_to(callstack); // Anonymous scope
              _run(ncommands_line,nposition,images,images_names,images,images_names,variables_sizes,0,0,0,false);
              callstack.remove();

            } else { // Not found -> Try generic image saver

              g_list.assign(selection.height());
              cimg_forY(selection,l) if (!gmic_check(images[selection(l)]))
                CImg<unsigned int>::vector(selection(l)).move_to(empty_indices);
              if (empty_indices && is_verbose) {
                selection2string(empty_indices>'y',images_names,1,eselec);
                warn(images,0,false,
                     "Command 'output': Image%s %s empty.",
                     eselec.data(),empty_indices.size()>1?"are":"is");
              }
              cimg_forY(selection,l)
                g_list[l].assign(images[selection[l]],true);
              if (g_list.size()==1)
                print(images,0,"Output image%s as %s file '%s' (1 image %dx%dx%dx%d).",
                      gmic_selection.data(),
                      uext.data(),_filename.data(),
                      g_list[0].width(),g_list[0].height(),
                      g_list[0].depth(),g_list[0].spectrum());
              else print(images,0,"Output image%s as %s file '%s'.",
                         gmic_selection.data(),uext.data(),_filename.data());

              if (*options)
                error(true,images,0,0,
                      "Command 'output': File '%s', format '%s' does not take any output options "
                      "(options '%s' specified).",
                      _filename.data(),ext,options.data());
              if (g_list.size()==1) g_list[0].save(filename); else g_list.save(filename);
            }
          }

          if (*cext) { // When output forced to 'ext' : copy final file to specified location
            try {
              CImg<unsigned char>::get_load_raw(filename_tmp).save_raw(_filename);
              std::remove(filename_tmp);
            } catch (...) { // Failed, maybe 'filename_tmp' consists of several numbered images
              bool save_failure = false;
              *message = 0;
              for (unsigned int i = 0; i!=~0U; ++i) {
                cimg::number_filename(filename_tmp,i,6,gmic_use_formula);
                cimg::number_filename(_filename,i,6,gmic_use_message);
                try { CImg<unsigned char>::get_load_raw(formula).save_raw(message); }
                catch (...) { i = ~0U - 1; if (!i) save_failure = true; }
              }
              if (save_failure)
                error(true,images,0,0,
                      "Command 'output': Invalid write of file '%s' from temporary file '%s'.",
                      _filename.data(),filename_tmp.data());
            }
          }
          g_list.assign();
          g_list_c.assign();
          if (is_stdout) std::fflush(stdout);
          is_change = false;
          ++position;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'p...'
        //-----------------------------
      gmic_commands_p :

        // Pass image from parent context.
        if (!std::strcmp("pass",command)) {
          gmic_substitute_args(false);
          if (cimg_sscanf(argument,"%d%c",&(err=2),&end)==1 && err>=-1 && err<=2) ++position;
          else err = 2;
          if (err<0)
            print(images,0,"Return ind%s of image%s from parent context.",
                  selection.height()==1?"ex":"ices",
                  gmic_selection.data());
          else {
            print(images,0,"Insert image%s from parent context %s%s.",
                  gmic_selection.data(),
                  err==0?"in non-shared state":
                  err==1?"in shared state":"using adaptive state",
                  selection.height()>1?"s":"");

            cimg_forY(selection,l) {
              CImg<T> &img = parent_images[selection[l]];
              const T *p = 0;
              std::memcpy(&p,&img._width,sizeof(void*));

              if (p && !img.data()) {
                // Parent image is in the current selection -> must search the current list
                bool found_image = false;
                cimglist_for(images,i) {
                  if (images[i].data()==p) { // Found it !
                    images.insert(images[i],~0U,err==1);
                    images_names[i].get_copymark().move_to(images_names);
                    found_image = true;
                    break;
                  }
                }
                if (!found_image) error(true,images,0,0,
                                        "Command 'pass': Unreferenced image [%d] from parent context "
                                        "(has been re-allocated in current context or reserved by another thread).",
                                        selection[l]);
              } else { // Parent image not in the current selection
                images.insert(img,~0U,(bool)err);
                images_names.insert(parent_images_names[selection[l]]);
              }
            }
          }
          selection.value_string().move_to(status);
          is_change = true;
          continue;
        }

        // Run multiple commands in parallel.
        if (!is_get && !std::strcmp("parallel",item)) {
          gmic_substitute_args(false);
          const char *_arg = argument, *_arg_text = gmic_argument_text_printed();
          bool wait_mode = true;
          if ((*_arg=='0' || *_arg=='1') && (_arg[1]==',' || !_arg[1])) {
            wait_mode = (bool)(*_arg - '0');
            uind = _arg[1]?2:1;
            _arg+=uind; _arg_text+=uind;
          }
          if (!*_arg) print(images,0,"Skip command 'parallel' (no commands specified)."); // No command specified
          else {
            CImgList<char> arguments = CImg<char>::string(_arg).get_split(CImg<char>::vector(','),0,false);
            CImg<_gmic_parallel<T> >(1,arguments.width()).move_to(gmic_threads);
            CImg<_gmic_parallel<T> > &_gmic_threads = gmic_threads.back();

#ifdef gmic_is_parallel
            print(images,0,"Execute %d command%s '%s' in parallel%s.",
                  arguments.width(),arguments.width()>1?"s":"",_arg_text,
                  wait_mode?" and wait thread termination immediately":
                  " and wait thread termination when current environment ends");
#else // #ifdef gmic_is_parallel
            print(images,0,"Execute %d commands '%s' (run sequentially, "
                  "parallel computing disabled).",
                  arguments.width(),_arg_text);
#endif // #ifdef gmic_is_parallel

            // Prepare thread structures.
            cimg_forY(_gmic_threads,l) {
              gmic &gmic_instance = _gmic_threads[l].gmic_instance;
              for (unsigned int i = 0; i<gmic_comslots; ++i) {
                gmic_instance.commands[i].assign(commands[i],true);
                gmic_instance.commands_names[i].assign(commands_names[i],true);
                gmic_instance.commands_has_arguments[i].assign(commands_has_arguments[i],true);
              }
              for (unsigned int i = 0; i<gmic_varslots; ++i) {
                if (i>=6*gmic_varslots/7) { // Share inter-thread global variables
                  gmic_instance.variables[i] = variables[i];
                  gmic_instance.variables_names[i] = variables_names[i];
                  gmic_instance.variables_lengths[i] = variables_lengths[i];
                } else {
                  if (i>=gmic_varslots/2) { // Make a copy of single-thread global variables
                    gmic_instance._variables[i].assign(_variables[i]);
                    gmic_instance._variables_names[i].assign(_variables_names[i]);
                    gmic_instance._variables_lengths[i].assign(_variables_lengths[i]);
                    _gmic_threads[l].variables_sizes[i] = variables_sizes[i];
                  } else _gmic_threads[l].variables_sizes[i] = 0;
                  gmic_instance.variables[i] = &gmic_instance._variables[i];
                  gmic_instance.variables_names[i] = &gmic_instance._variables_names[i];
                  gmic_instance.variables_lengths[i] = &gmic_instance._variables_lengths[i];
                }
              }
              gmic_instance.callstack.assign(callstack);
              gmic_instance.commands_files.assign(commands_files,true);
              gmic_use_title;
              cimg_snprintf(title,_title.width(),"*thread%d",l);
              CImg<char>::string(title).move_to(gmic_instance.callstack);
              gmic_instance.light3d.assign(light3d);
              gmic_instance.status.assign(status);
              gmic_instance.debug_filename = debug_filename;
              gmic_instance.debug_line = debug_line;
              gmic_instance.light3d_x = light3d_x;
              gmic_instance.light3d_y = light3d_y;
              gmic_instance.light3d_z = light3d_z;
              gmic_instance._progress = 0;
              gmic_instance.progress = &gmic_instance._progress;
              gmic_instance.is_change = is_change;
              gmic_instance.is_debug = is_debug;
              gmic_instance.is_start = false;
              gmic_instance.is_quit = false;
              gmic_instance.is_return = false;
              gmic_instance.verbosity = verbosity;
              gmic_instance._is_abort = _is_abort;
              gmic_instance.is_abort = is_abort;
              gmic_instance.is_abort_thread = false;
              gmic_instance.nb_carriages_default = nb_carriages_default;
              gmic_instance.nb_carriages_stdout = nb_carriages_stdout;
              gmic_instance.reference_time = reference_time;
              _gmic_threads[l].images = &images;
              _gmic_threads[l].images_names = &images_names;
              _gmic_threads[l].parent_images = &parent_images;
              _gmic_threads[l].parent_images_names = &parent_images_names;
              _gmic_threads[l].gmic_threads = &_gmic_threads;
              _gmic_threads[l].command_selection = command_selection;
              _gmic_threads[l].is_thread_running = true;

              // Substitute special characters codes appearing outside strings.
              arguments[l].resize(1,arguments[l].height() + 1,1,1,0);
              bool is_dquoted = false;
              for (char *s = arguments[l].data(); *s; ++s) {
                if (*s=='\"') is_dquoted = !is_dquoted;
                else if (!is_dquoted) _strreplace_fw(*s);
              }
              gmic_instance.commands_line_to_CImgList(arguments[l].data()).
                move_to(_gmic_threads[l].commands_line);
            }

            // Run threads.
            cimg_forY(_gmic_threads,l) {
#ifdef gmic_is_parallel
#ifdef PTHREAD_CANCEL_ENABLE

#if defined(__MACOSX__) || defined(__APPLE__)
              const cimg_uint64 stacksize = (cimg_uint64)8*1024*1024;
              pthread_attr_t thread_attr;
              if (!pthread_attr_init(&thread_attr) && !pthread_attr_setstacksize(&thread_attr,stacksize))
                // Reserve enough stack size for the new thread.
                pthread_create(&_gmic_threads[l].thread_id,&thread_attr,gmic_parallel<T>,(void*)&_gmic_threads[l]);
              else
#endif // #if defined(__MACOSX__) || defined(__APPLE__)
                pthread_create(&_gmic_threads[l].thread_id,0,gmic_parallel<T>,(void*)&_gmic_threads[l]);

#elif cimg_OS==2 // #ifdef PTHREAD_CANCEL_ENABLE
              _gmic_threads[l].thread_id = CreateThread(0,0,(LPTHREAD_START_ROUTINE)gmic_parallel<T>,
                                                        (void*)&_gmic_threads[l],0,0);
#endif // #ifdef PTHREAD_CANCEL_ENABLE
#else // #ifdef gmic_is_parallel
              gmic_parallel<T>((void*)&_gmic_threads[l]);
#endif // #ifdef gmic_is_parallel
            }

            // Wait threads if immediate waiting mode selected.
            if (wait_mode) {
              wait_threads((void*)&_gmic_threads,false,(T)0);

              // Get 'released' state of the image list.
              cimg_forY(_gmic_threads,l) is_change|=_gmic_threads[l].gmic_instance.is_change;

              // Get status modified by first thread.
              _gmic_threads[0].gmic_instance.status.move_to(status);

              // Check for possible exceptions thrown by threads.
              cimg_forY(_gmic_threads,l) if (_gmic_threads[l].exception._message)
                error(false,images,0,_gmic_threads[l].exception.command(),"%s",_gmic_threads[l].exception.what());
              gmic_threads.remove();
            }
          }
          ++position;
          continue;
        }

        // Permute axes.
        if (!std::strcmp("permute",command)) {
          gmic_substitute_args(false);
          print(images,0,"Permute axes of image%s, with permutation '%s'.",
                gmic_selection.data(),gmic_argument_text_printed());
          cimg_forY(selection,l) gmic_apply(permute_axes(argument),false);
          is_change = true;
          ++position;
          continue;
        }

        // Set progress index.
        if (!is_get && !std::strcmp("progress",item)) {
          gmic_substitute_args(false);
          value = -1;
          if (cimg_sscanf(argument,"%lf%c",&value,&end)!=1) {
            name.assign(argument,(unsigned int)std::strlen(argument) + 1);
            CImg<T> &img = images.size()?images.back():CImg<T>::empty();
            strreplace_fw(name);
            try { value = img.eval(name,0,0,0,0,&images); }
            catch (CImgException &e) {
              const char *const e_ptr = std::strstr(e.what(),": ");
              error(true,images,0,"progress",
                    "Command 'progress': Invalid argument '%s': %s",
                    cimg::strellipsize(name,64,false),e_ptr?e_ptr + 2:e.what());
            }
          }
          if (value<0) value = -1; else if (value>100) value = 100;
          if (value>=0)
            print(images,0,"Set progress index to %g%%.",
                  value);
          else
            print(images,0,"Disable progress index.");
          *progress = (float)value;
          ++position;
          continue;
        }

        // Print.
        if (!std::strcmp("print",command)) {
          const int _verbosity = ++verbosity;
          std::FILE *_file = 0;
          if (is_get) { _file = cimg::output(); verbosity = 1; cimg::output(stdout); }
          print_images(images,images_names,selection);
          if (is_get) { verbosity = _verbosity; cimg::output(_file); }
          --verbosity;
          is_change = false;
          continue;
        }

        // Power.
        gmic_arithmetic_command("pow",
                                pow,
                                "Compute image%s to the power of %g%s",
                                gmic_selection.data(),value,ssep,Tfloat,
                                pow,
                                "Compute image%s to the power of image [%d]",
                                gmic_selection.data(),ind[0],
                                pow,
                                "Compute image%s to the power of expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute sequential power of image%s");

        // Draw point.
        if (!std::strcmp("point",command)) {
          gmic_substitute_args(false);
          double x = 0, y = 0, z = 0;
          sepx = sepy = sepz = *argx = *argy = *argz = *color = 0;
          opacity = 1;
          if ((cimg_sscanf(argument,"%255[0-9.eE%+-]%c",
                           gmic_use_argx,&end)==1 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           argx,gmic_use_argy,&end)==2 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           argx,argy,gmic_use_argz,&end)==3 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%f%c",
                           argx,argy,argz,&opacity,&end)==4 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%f,"
                           "%4095[0-9.eEinfa,+-]%c",
                           argx,argy,argz,&opacity,gmic_use_color,&end)==5) &&
              (cimg_sscanf(argx,"%lf%c",&x,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy ||
               cimg_sscanf(argy,"%lf%c",&y,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
              (!*argz ||
               cimg_sscanf(argz,"%lf%c",&z,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&z,&sepz,&end)==2 && sepz=='%'))) {
            print(images,0,
                  "Draw point (%g%s,%g%s,%g%s) on image%s, with opacity %g and color (%s).",
                  x,sepx=='%'?"%":"",
                  y,sepy=='%'?"%":"",
                  z,sepz=='%'?"%":"",
                  gmic_selection.data(),
                  opacity,
                  *color?color:"default");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              g_img.assign(img.spectrum(),1,1,1,(T)0).fill_from_values(color,true);
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.width() - 1)/100:x),
                ny = (int)cimg::round(sepy=='%'?y*(img.height() - 1)/100:y),
                nz = (int)cimg::round(sepz=='%'?z*(img.depth() - 1)/100:z);
              gmic_apply(draw_point(nx,ny,nz,g_img.data(),opacity),true);
            }
          } else arg_error("point");
          g_img.assign();
          is_change = true;
          ++position;
          continue;
        }

        // Draw polygon.
        if (!std::strcmp("polygon",command)) {
          gmic_substitute_args(false);
          name.assign(256);
          double N = 0, x0 = 0, y0 = 0;
          sep1 = sepx = sepy = *name = *color = 0;
          pattern = ~0U; opacity = 1;
          if (cimg_sscanf(argument,"%lf%c",
                          &N,&end)==2 && N>=1) {
            N = cimg::round(N);
            const char
              *nargument = argument + cimg_snprintf(name,name.width(),"%u",
                                                    (unsigned int)N) + 1,
              *const eargument = argument + std::strlen(argument);
            vertices.assign((unsigned int)N,2,1,1,0);
            CImg<bool> percents((unsigned int)N,2,1,1,0);
            for (unsigned int n = 0; n<(unsigned int)N; ++n) if (nargument<eargument) {
                sepx = sepy = 0;
                if (cimg_sscanf(nargument,"%255[0-9.eE%+-],%255[0-9.eE%+-]",
                                gmic_use_argx,gmic_use_argy)==2 &&
                    (cimg_sscanf(argx,"%lf%c",&x0,&end)==1 ||
                     (cimg_sscanf(argx,"%lf%c%c",&x0,&sepx,&end)==2 && sepx=='%')) &&
                    (cimg_sscanf(argy,"%lf%c",&y0,&end)==1 ||
                     (cimg_sscanf(argy,"%lf%c%c",&y0,&sepy,&end)==2 && sepy=='%'))) {
                  vertices(n,0U) = (float)x0; percents(n,0U) = (sepx=='%');
                  vertices(n,1U) = (float)y0; percents(n,1U) = (sepy=='%');
                  nargument+=std::strlen(argx) + std::strlen(argy) + 2;
                } else arg_error("polygon");
              } else arg_error("polygon");
            if (nargument<eargument &&
                cimg_sscanf(nargument,"%4095[0-9.eEinfa+-]",gmic_use_color)==1 &&
                cimg_sscanf(color,"%f",&opacity)==1) {
              nargument+=std::strlen(color) + 1;
              *color = 0;
            }
            if (nargument<eargument &&
                cimg_sscanf(nargument,"0%c%4095[0-9a-fA-F]",&sep1,gmic_use_color)==2 && sep1=='x' &&
                cimg_sscanf(color,"%x%c",&pattern,&end)==1) {
              nargument+=std::strlen(color) + 3;
              *color = 0;
            }
            const char *const p_color = nargument<eargument?nargument:&(end=0);
            if (sep1=='x')
              print(images,0,"Draw %g-vertices outlined polygon on image%s, with opacity %g, "
                    "pattern 0x%x and color (%s).",
                    N,
                    gmic_selection.data(),
                    opacity,pattern,
                    *p_color?p_color:"default");
            else
              print(images,0,"Draw %g-vertices filled polygon on image%s, with opacity %g "
                    "and color (%s).",
                    N,
                    gmic_selection.data(),
                    opacity,
                    *p_color?p_color:"default");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              CImg<int> coords(vertices.width(),2,1,1,0);
              cimg_forX(coords,p) {
                if (percents(p,0))
                  coords(p,0) = (int)cimg::round(vertices(p,0)*(img.width() - 1)/100);
                else coords(p,0) = (int)cimg::round(vertices(p,0));
                if (percents(p,1))
                  coords(p,1) = (int)cimg::round(vertices(p,1)*(img.height() - 1)/100);
                else coords(p,1) = (int)cimg::round(vertices(p,1));
              }
              g_img.assign(img.spectrum(),1,1,1,(T)0).fill_from_values(p_color,true);
              if (sep1=='x') { gmic_apply(draw_polygon(coords,g_img.data(),opacity,pattern),true); }
              else gmic_apply(draw_polygon(coords,g_img.data(),opacity),true);
            }
          } else arg_error("polygon");
          vertices.assign();
          g_img.assign();
          is_change = true;
          ++position;
          continue;
        }

        // Draw plasma fractal.
        if (!std::strcmp("plasma",command)) {
          gmic_substitute_args(false);
          float alpha = 1, beta = 1, scale = 8;
          if ((cimg_sscanf(argument,"%f%c",
                           &alpha,&end)==1 ||
               cimg_sscanf(argument,"%f,%f%c",
                           &alpha,&beta,&end)==2 ||
               cimg_sscanf(argument,"%f,%f,%f%c",
                           &alpha,&beta,&scale,&end)==3) &&
              scale>=0) ++position;
          else { alpha = beta = 1; scale = 8; }
          const unsigned int _scale = (unsigned int)cimg::round(scale);
          print(images,0,"Draw plasma fractal on image%s, with alpha %g, beta %g and scale %u.",
                gmic_selection.data(),
                alpha,
                beta,
                _scale);
          cimg_forY(selection,l) gmic_apply(draw_plasma(alpha,beta,_scale),false);
          is_change = true;
          continue;
        }

        // Display as a graph plot.
        if (!is_get && !std::strcmp("plot",command)) {
          gmic_substitute_args(false);
          double ymin = 0, ymax = 0, xmin = 0, xmax = 0, resolution = 65536;
          unsigned int plot_type = 1, vertex_type = 1;
          *formula = sep = 0;
          exit_on_anykey = 0;
          if (((cimg_sscanf(argument,"'%1023[^']%c%c",
                            gmic_use_formula,&sep,&end)==2 && sep=='\'') ||
               cimg_sscanf(argument,"'%1023[^']',%lf%c",
                           formula,&resolution,&end)==2 ||
               cimg_sscanf(argument,"'%1023[^']',%lf,%u%c",
                           formula,&resolution,&plot_type,&end)==3 ||
               cimg_sscanf(argument,"'%1023[^']',%lf,%u,%u%c",
                           formula,&resolution,&plot_type,&vertex_type,&end)==4 ||
               cimg_sscanf(argument,"'%1023[^']',%lf,%u,%u,%lf,%lf%c",
                           formula,&resolution,&plot_type,&vertex_type,&xmin,&xmax,&end)==6 ||
               cimg_sscanf(argument,"'%1023[^']',%lf,%u,%u,%lf,%lf,%lf,%lf%c",
                           formula,&resolution,&plot_type,&vertex_type,
                           &xmin,&xmax,&ymin,&ymax,&end)==8 ||
               cimg_sscanf(argument,"'%1023[^']',%lf,%u,%u,%lf,%lf,%lf,%lf,%u%c",
                           formula,&resolution,&plot_type,&vertex_type,
                           &xmin,&xmax,&ymin,&ymax,&exit_on_anykey,&end)==9) &&
              resolution>0 && plot_type<=3 && vertex_type<=7 && exit_on_anykey<=1) {
            resolution = cimg::round(resolution);
            strreplace_fw(formula);
            if (xmin==0 && xmax==0) { xmin = -4; xmax = 4; }
            if (!plot_type && !vertex_type) plot_type = 1;
            if (resolution<1) resolution = 65536;

            gmic_use_argx;
            cimg_snprintf(argx,_argx.width(),"x = lerp(%g,%g,x/%d);",
                          xmin,xmax,(unsigned int)(resolution>1?resolution - 1:0));
            const CImg<char> n_formula = (CImg<char>::string(argx,false,true),
                                          CImg<char>::string(formula,true,true))>'x';
            boundary = 1U;
            try { // Determine vector dimension of specified formula
              typename CImg<T>::_cimg_math_parser mp(n_formula.data() + (*n_formula=='>' || *n_formula=='<' ||
                                                                         *n_formula=='*' || *n_formula==':'),
                                                     "plot",CImg<T>::const_empty(),0,&images);
              boundary = std::max(1U,mp.result_dim);
            } catch (...) { is_cond = false; }
            CImg<T> values((int)resolution,1,1,boundary,0);
            values.fill(n_formula,false,true,&images);

            gmic_use_title;
            cimg_snprintf(title,_title.width(),"[Plot of '%s']",formula);
            CImg<char>::string(title).move_to(g_list_c);
            display_plots(CImgList<T>(values,true),g_list_c,CImg<unsigned int>::vector(0),
                          plot_type,vertex_type,xmin,xmax,ymin,ymax,exit_on_anykey);
            g_list_c.assign();
            ++position;
          } else {
            plot_type = 1; vertex_type = 0; ymin = ymax = xmin = xmax = 0;
            if ((cimg_sscanf(argument,"%u%c",
                             &plot_type,&end)==1 ||
                 cimg_sscanf(argument,"%u,%u%c",
                             &plot_type,&vertex_type,&end)==2 ||
                 cimg_sscanf(argument,"%u,%u,%lf,%lf%c",
                             &plot_type,&vertex_type,&xmin,&xmax,&end)==4 ||
                 cimg_sscanf(argument,"%u,%u,%lf,%lf,%lf,%lf%c",
                             &plot_type,&vertex_type,&xmin,&xmax,&ymin,&ymax,&end)==6 ||
                 cimg_sscanf(argument,"%u,%u,%lf,%lf,%lf,%lf,%u%c",
                             &plot_type,&vertex_type,&xmin,&xmax,&ymin,&ymax,&exit_on_anykey,&end)==7) &&
                plot_type<=3 && vertex_type<=7 && exit_on_anykey<=1) ++position;
            if (!plot_type && !vertex_type) plot_type = 1;
            display_plots(images,images_names,selection,plot_type,vertex_type,
                          xmin,xmax,ymin,ymax,exit_on_anykey);
          }
          is_change = false;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'q...'
        //-----------------------------
      gmic_commands_q :

        // Quit.
        if (!is_get && !std::strcmp("quit",item)) {
          print(images,0,"Quit G'MIC interpreter.");
          dowhiles.assign(nb_dowhiles = 0U);
          fordones.assign(nb_fordones = 0U);
          foreachdones.assign(nb_foreachdones = 0U);
          repeatdones.assign(nb_repeatdones = 0U);
          position = commands_line.size();
          is_change = false;
          is_quit = true;
          *is_abort = true;
          break;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'r...'
        //-----------------------------
      gmic_commands_r :

        // Remove images.
        if (!std::strcmp("remove",command)) {
          print(images,0,"Remove image%s",
                gmic_selection.data());
          if (is_get) {
            g_list.assign(images);
            g_list_c.assign(images_names);
          }
          remove_images(images,images_names,selection,0,selection.height() - 1);
          if (is_get) {
            cimglist_for(images,l) images_names[l].get_copymark().move_to(images_names[l]);
            g_list.move_to(images,0);
            g_list_c.move_to(images_names,0);
          }
          if (is_verbose) {
            cimg::mutex(29);
            std::fprintf(cimg::output()," (%u image%s left).",
                         images.size(),images.size()==1?"":"s");
            std::fflush(cimg::output());
            cimg::mutex(29,0);
          }
          g_list.assign();
          g_list_c.assign();
          is_change = true;
          continue;
        }

        // Repeat.
        if (!is_get && !std::strcmp("repeat",item)) {
          gmic_substitute_args(false);
          if (cimg_sscanf(argument,"%lf%c",&value,&end)!=1) {
            name.assign(argument,(unsigned int)std::strlen(argument) + 1);
            strreplace_fw(name);
            CImg<T> &img = images.size()?images.back():CImg<T>::empty();
            try { value = img.eval(name,0,0,0,0,&images); }
            catch (CImgException &e) {
              const char *const e_ptr = std::strstr(e.what(),": ");
              error(true,images,0,"repeat",
                    "Command 'repeat': Invalid argument '%s'; %s",
                    cimg::strellipsize(name,64,false),e_ptr?e_ptr + 2:e.what());
            }
          }
          const unsigned int nb = value<=0?0U:
            cimg::type<double>::is_inf(value)?~0U:(unsigned int)cimg::round(value);
          ++position;

          if (!nb) {
            if (is_very_verbose) print(images,0,"Skip 'repeat...done' block (0 iterations).");
            int nb_levels = 0;
            for (nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
              it = commands_line[position];
              if (*it==1)
                is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
              else {
                _is_get = *it=='+';
                it+=(_is_get || *it=='-');
                gmic_if_flr ++nb_levels; gmic_elif_flr --nb_levels;
              }
            }
            if (nb_levels)
              error(true,images,0,0,
                    "Command 'repeat': Missing associated 'done' command.");
          } else {
            if (is_debug_info && debug_line!=~0U) {
              gmic_use_argx;
              cimg_snprintf(argx,_argx.width(),"*repeat#%u",debug_line);
              CImg<char>::string(argx).move_to(callstack);
            } else CImg<char>::string("*repeat").move_to(callstack);
            if (is_very_verbose) print(images,0,"Start 'repeat...done' block (%u iteration%s).",nb,nb>1?"s":"");
            if (nb_repeatdones>=repeatdones._height) repeatdones.resize(4,std::max(2*repeatdones._height,8U),1,1,0);
            unsigned int *const rd = repeatdones.data(0,nb_repeatdones++);
            rd[0] = position_item;
            rd[1] = 0;
            rd[2] = nb;
            rd[3] = debug_line;
          }
          is_lbrace_command = true;
          continue;
        }

        // Resize.
        if (!std::strcmp("resize",command)) {
          gmic_substitute_args(true);
          double valx = 100, valy = 100, valz = 100, valc = 100;
          float cx = 0, cy = 0, cz = 0, cc = 0;
          CImg<char> indicesy(256), indicesz(256), indicesc(256);
          CImg<unsigned int> indx, indy, indz, indc;
          *indices = *indicesy = *indicesz = *indicesc = *argx = *argy = *argz = *argc = sep = 0;
          sepx = sepy = sepz = sepc = '%';
          int iinterpolation = 1;
          boundary = 0;
          ind.assign();
          if ((cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-]%c",
                           gmic_use_argx,&end)==1 ||
               cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-]%c",
                           argx,gmic_use_argy,&end)==2 ||
               cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],"
                           "%255[][a-zA-Z0-9_.eE%+-]%c",
                           argx,argy,gmic_use_argz,&end)==3 ||
               cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],"
                           "%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-]%c",
                           argx,argy,argz,gmic_use_argc,&end)==4 ||
               cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],"
                           "%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],%d%c",
                           argx,argy,argz,argc,&iinterpolation,&end)==5 ||
               cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],"
                           "%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],%d,%u%c",
                           argx,argy,argz,argc,&iinterpolation,&boundary,&end)==6 ||
               cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],"
                           "%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],%d,%u,%f%c",
                           argx,argy,argz,argc,&iinterpolation,&boundary,&cx,&end)==7 ||
               cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],"
                           "%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],%d,%u,%f,"
                           "%f%c",
                           argx,argy,argz,argc,&iinterpolation,&boundary,&cx,&cy,&end)==8||
               cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],"
                           "%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],%d,%u,%f,"
                           "%f,%f%c",
                           argx,argy,argz,argc,&iinterpolation,&boundary,
                           &cx,&cy,&cz,&end)==9 ||
               cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],"
                           "%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],%d,%u,%f,"
                           "%f,%f,%f%c",
                           argx,argy,argz,argc,&iinterpolation,&boundary,
                           &cx,&cy,&cz,&cc,&end)==10) &&
              ((cimg_sscanf(argx,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sepx,&end)==2 &&
                sepx==']' &&
                (indx=selection2cimg(indices,images.size(),images_names,"resize")).height()==1) ||
               (sepx=0,cimg_sscanf(argx,"%lf%c",&valx,&sepx)==1 && valx>=1) ||
               (cimg_sscanf(argx,"%lf%c%c",&valx,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy ||
               (cimg_sscanf(argy,"[%255[a-zA-Z0-9_.%+-]%c%c",indicesy.data(),&sepy,
                            &end)==2 &&
                sepy==']' &&
                (indy=selection2cimg(indicesy,images.size(),images_names,"resize")).height()==1) ||
               (sepy=0,cimg_sscanf(argy,"%lf%c",&valy,&sepy)==1 && valy>=1) ||
               (cimg_sscanf(argy,"%lf%c%c",&valy,&sepy,&end)==2 && sepy=='%')) &&
              (!*argz ||
               (cimg_sscanf(argz,"[%255[a-zA-Z0-9_.%+-]%c%c",indicesz.data(),&sepz,
                            &end)==2 &&
                sepz==']' &&
                (indz=selection2cimg(indicesz,images.size(),images_names,"resize")).height()==1) ||
               (sepz=0,cimg_sscanf(argz,"%lf%c",&valz,&sepz)==1 && valz>=1) ||
               (cimg_sscanf(argz,"%lf%c%c",&valz,&sepz,&end)==2 && sepz=='%')) &&
              (!*argc ||
               (cimg_sscanf(argc,"[%255[a-zA-Z0-9_.%+-]%c%c",indicesc.data(),&sepc,
                            &end)==2 &&
                sepc==']' &&
                (indc=selection2cimg(indicesc,images.size(),images_names,"resize")).height()==1) ||
               (sepc=0,cimg_sscanf(argc,"%lf%c",&valc,&sepc)==1 && valc>=1) ||
               (cimg_sscanf(argc,"%lf%c%c",&valc,&sepc,&end)==2 && sepc=='%')) &&
              valx>0 && valy>0 && valz>0 && valc>0 &&
              iinterpolation>=-1 && iinterpolation<=6 && boundary<=3 &&
              cx>=0 && cx<=1 && cy>=0 && cy<=1 && cz>=0 && cz<=1 && cc>=0 && cc<=1) {
            if (indx) { valx = (float)images[*indx].width(); sepx = 0; }
            if (indy) { valy = (float)images[*indy].height(); sepy = 0; }
            if (indz) { valz = (float)images[*indz].depth(); sepz = 0; }
            if (indc) { valc = (float)images[*indc].spectrum(); sepc = 0; }
            print(images,0,
                  "Resize image%s to %g%s%g%s%g%s%g%s, with %s interpolation, "
                  "%s boundary conditions and alignment (%g,%g,%g,%g).",
                  gmic_selection.data(),
                  valx,sepx=='%'?"%x":"x",
                  valy,sepy=='%'?"%x":"x",
                  valz,sepz=='%'?"%x":"x",
                  valc,sepc=='%'?"% ":"",
                  iinterpolation<=0?"no":iinterpolation==1?"nearest-neighbor":
                  iinterpolation==2?"moving average":iinterpolation==3?"linear":
                  iinterpolation==4?"grid":iinterpolation==5?"cubic":"lanczos",
                  boundary<=0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror",
                  cx,cy,cz,cc);
            cimg_forY(selection,l) {
              CImg<T>& img = images[selection[l]];
              const int
                _nvalx = (int)cimg::round(sepx=='%'?valx*img.width()/100:valx),
                _nvaly = (int)cimg::round(sepy=='%'?valy*img.height()/100:valy),
                _nvalz = (int)cimg::round(sepz=='%'?valz*img.depth()/100:valz),
                _nvalc = (int)cimg::round(sepc=='%'?valc*img.spectrum()/100:valc),
                nvalx = _nvalx?_nvalx:1,
                nvaly = _nvaly?_nvaly:1,
                nvalz = _nvalz?_nvalz:1,
                nvalc = _nvalc?_nvalc:1;
              gmic_apply(resize(nvalx,nvaly,nvalz,nvalc,iinterpolation,boundary,cx,cy,cz,cc),false);
            }
          } else arg_error("resize");
          is_change = true;
          ++position;
          continue;
        }

        // Reverse positions.
        if (!std::strcmp("reverse",command)) {
          print(images,0,"Reverse positions of image%s.",
                gmic_selection.data());
          if (is_get) {
            pattern = images.size();
            images.insert(selection.height());
            images_names.insert(selection.height());
            cimg_forY(selection,l) {
              uind = selection[selection.height() - 1 - l];
              images[pattern + l].assign(images[uind],false);
              images_names[uind].get_copymark().move_to(images_names[pattern + l]);
            }
          } else for (int l = 0; l<selection.height()/2; ++l) {
              const unsigned int
                i0 = selection[l],
                i1 = selection[selection.height() - 1 - l];
              images[i0].swap(images[i1]);
              images_names[i0].swap(images_names[i1]);
            }
          is_change = true;
          continue;
        }

        // Return.
        if (!is_get && !std::strcmp("return",item)) {
          if (is_very_verbose) print(images,0,"Return.");
          position = commands_line.size();
          is_return = true;
          break;
        }

        // Rotate.
        if (!std::strcmp("rotate",command)) {
          gmic_substitute_args(false);
          float angle = 0, u = 0, v = 0, w = 0, cx = 0, cy = 0, cz = 0;
          char sep2 = sep1 = sep0 = *argx = *argy = *argz = 0;
          interpolation = 1;
          boundary = 0;
          if ((cimg_sscanf(argument,"%f%c",
                           &angle,&end)==1 ||
               cimg_sscanf(argument,"%f,%u%c",
                           &angle,&interpolation,&end)==2 ||
               cimg_sscanf(argument,"%f,%u,%u%c",
                           &angle,&interpolation,&boundary,&end)==3 ||
               cimg_sscanf(argument,"%f,%u,%u,%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           &angle,&interpolation,&boundary,gmic_use_argx,gmic_use_argy,&end)==5) &&
              (!*argx ||
               cimg_sscanf(argx,"%f%c",&cx,&end)==1 ||
               (cimg_sscanf(argx,"%f%c%c",&cx,&sep0,&end)==2 && sep0=='%')) &&
              (!*argy ||
               cimg_sscanf(argy,"%f%c",&cy,&end)==1 ||
               (cimg_sscanf(argy,"%f%c%c",&cy,&sep1,&end)==2 && sep1=='%')) &&
              interpolation<=2 && boundary<=3) { // 2D rotation
            if (*argx) {
              print(images,0,"Rotate image%s of %g deg., with %s interpolation, %s boundary conditions "
                    "and center at (%g%s,%g%s).",
                    gmic_selection.data(),angle,
                    interpolation==0?"nearest-neighbor":interpolation==1?"linear":"cubic",
                    boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror",
                    cx,sep0=='%'?"%":"",cy,sep1=='%'?"%":"");
              cimg_forY(selection,l) {
                CImg<T> &img = images[selection[l]];
                const float
                  ncx = sep0=='%'?cx*(img.width() - 1)/100:cx,
                  ncy = sep1=='%'?cy*(img.height() - 1)/100:cy;
                gmic_apply(rotate(angle,ncx,ncy,interpolation,boundary),false);
              }
            } else {
              print(images,0,"Rotate image%s of %g deg., with %s interpolation and %s boundary conditions.",
                    gmic_selection.data(),angle,
                    interpolation==0?"nearest-neighbor":interpolation==1?"linear":"cubic",
                    boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror");
              cimg_forY(selection,l) gmic_apply(rotate(angle,interpolation,boundary),false);
            }
          } else if ((cimg_sscanf(argument,"%f,%f,%f,%f,%u,%u,%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                                  &u,&v,&w,&angle,&interpolation,&boundary,
                                  &(*gmic_use_argx=0),&(*gmic_use_argy=0),&(*gmic_use_argz=0),&end)==9 ||
                      cimg_sscanf(argument,"%f,%f,%f,%f,%u,%u%c",
                                  &u,&v,&w,&angle,&interpolation,&boundary,&end)==6) &&
                     (!*argx ||
                      cimg_sscanf(argx,"%f%c",&cx,&end)==1 ||
                      (cimg_sscanf(argx,"%f%c%c",&cx,&sep0,&end)==2 && sep0=='%')) &&
                     (!*argy ||
                      cimg_sscanf(argy,"%f%c",&cy,&end)==1 ||
                      (cimg_sscanf(argy,"%f%c%c",&cy,&sep1,&end)==2 && sep1=='%')) &&
                     (!*argz ||
                      cimg_sscanf(argz,"%f%c",&cz,&end)==1 ||
                      (cimg_sscanf(argz,"%f%c%c",&cz,&sep2,&end)==2 && sep2=='%')) &&
                     interpolation<=2 && boundary<=3) { // 3D rotation
            if (*argx) {
              print(images,0,"Rotate image%s around axis (%g,%g,%g) with angle %g deg., %s interpolation, "
                    "%s boundary conditions and center at (%g%s,%g%s,%g%s).",
                    gmic_selection.data(),u,v,w,angle,
                    interpolation==0?"nearest-neighbor":interpolation==1?"linear":"cubic",
                    boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror",
                    cx,sep0=='%'?"%":"",cy,sep1=='%'?"%":"",cz,sep2=='%'?"%":"");
              cimg_forY(selection,l) {
                CImg<T> &img = images[selection[l]];
                const float
                  ncx = sep0=='%'?cx*(img.width() - 1)/100:cx,
                  ncy = sep1=='%'?cy*(img.height() - 1)/100:cy,
                  ncz = sep2=='%'?cy*(img.depth() - 1)/100:cz;
                gmic_apply(rotate(u,v,w,angle,ncx,ncy,ncz,interpolation,boundary),false);
              }
            } else {
              print(images,0,"Rotate image%s around axis (%g,%g,%g) with angle %g deg., %s interpolation "
                    "and %s boundary conditions.",
                    gmic_selection.data(),u,v,w,angle,
                    interpolation==0?"nearest-neighbor":interpolation==1?"linear":"cubic",
                    boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror");
              cimg_forY(selection,l) gmic_apply(rotate(u,v,w,angle,interpolation,boundary),false);
            }
          } else arg_error("rotate");
          is_change = true;
          ++position;
          continue;
        }

        // Round.
        if (!std::strcmp("round",command)) {
          gmic_substitute_args(false);
          int rounding_type = 0;
          value = 1;
          if ((cimg_sscanf(argument,"%lf%c",
                           &value,&end)==1 ||
               cimg_sscanf(argument,"%lf,%d%c",
                           &value,&rounding_type,&end)==2) &&
              value>=0 && rounding_type>=-1 && rounding_type<=1) ++position;
          else { value = 1; rounding_type = 0; }
          print(images,0,"Round values of image%s by %g and %s rounding.",
                gmic_selection.data(),
                value,
                rounding_type<0?"backward":rounding_type>0?"forward":"nearest");
          cimg_forY(selection,l) gmic_apply(round(value,rounding_type),true);
          is_change = true;
          continue;
        }

        // Fill with random values.
        if (!std::strcmp("rand",command)) {
          gmic_substitute_args(true);
          ind0.assign(); ind1.assign();
          sep0 = sep1 = *argx = *argy = *indices = 0;
          value0 = value1 = 0;
          if (cimg_sscanf(argument,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-]%c",
                          gmic_use_argx,gmic_use_argy,&end)==2 &&
              ((cimg_sscanf(argx,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep0,&end)==2 &&
                sep0==']' &&
                (ind0=selection2cimg(indices,images.size(),images_names,"rand")).height()==1) ||
               (cimg_sscanf(argx,"%lf%c%c",&value0,&sep0,&end)==2 && sep0=='%') ||
               cimg_sscanf(argx,"%lf%c",&value0,&end)==1) &&
              ((cimg_sscanf(argy,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_formula,&sep1,&end)==2 &&
                sep1==']' &&
                (ind1=selection2cimg(formula,images.size(),images_names,"rand")).height()==1) ||
               (cimg_sscanf(argy,"%lf%c%c",&value1,&sep1,&end)==2 && sep1=='%') ||
               cimg_sscanf(argy,"%lf%c",&value1,&end)==1)) {
            if (ind0) { value0 = images[*ind0].min(); sep0 = 0; }
            if (ind1) { value1 = images[*ind1].max(); sep1 = 0; }
            print(images,0,"Fill image%s with random values, in range [%g%s,%g%s].",
                  gmic_selection.data(),
                  value0,sep0=='%'?"%":"",
                  value1,sep1=='%'?"%":"");
            cimg_forY(selection,l) {
              CImg<T>& img = gmic_check(images[selection[l]]);
              nvalue0 = value0; nvalue1 = value1;
              vmin = vmax = 0;
              if (sep0=='%' || sep1=='%') {
                if (img) vmax = (double)img.max_min(vmin);
                if (sep0=='%') nvalue0 = vmin + (vmax - vmin)*value0/100;
                if (sep1=='%') nvalue1 = vmin + (vmax - vmin)*value1/100;
              }
              gmic_apply(rand((T)nvalue0,(T)nvalue1),true);
            }
          } else if (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep0,&end)==2 &&
                     sep0==']' &&
                     (ind0=selection2cimg(indices,images.size(),images_names,"rand")).height()==1) {
            if (images[*ind0]) value1 = (double)images[*ind0].max_min(value0);
            print(images,0,"Fill image%s with random values, in range [%g,%g] from image [%d].",
                  gmic_selection.data(),
                  value0,
                  value1,
                  *ind0);
            cimg_forY(selection,l) gmic_apply(rand((T)value0,(T)value1),true);
          } else arg_error("rand");
          is_change = true;
          ++position;
          continue;
        }

        // Rotate 3D object.
        if (!std::strcmp("rotate3d",command)) {
          gmic_substitute_args(false);
          float u = 0, v = 0, w = 1, angle = 0;
          if (cimg_sscanf(argument,"%f,%f,%f,%f%c",
                          &u,&v,&w,&angle,&end)==4) {
            print(images,0,"Rotate 3D object%s around axis (%g,%g,%g), with angle %g deg.",
                  gmic_selection.data(),
                  u,v,w,
                  angle);
            CImg<float>::rotation_matrix(u,v,w,angle).move_to(vertices);
            cimg_forY(selection,l) {
              uind = selection[l];
              CImg<T>& img = images[uind];
              try { gmic_apply(rotate_CImg3d(vertices),true); }
              catch (CImgException&) {
                if (!img.is_CImg3d(true,&(*gmic_use_message=0)))
                  error(true,images,0,0,
                        "Command 'rotate3d': Invalid 3D object [%d], "
                        "in image%s (%s).",
                        uind,gmic_selection_err.data(),message);
                else throw;
              }
            }
          } else arg_error("rotate3d");
          vertices.assign();
          is_change = true;
          ++position;
          continue;
        }

        // Bitwise left rotation.
        gmic_arithmetic_command("rol",
                                rol,
                                "Compute bitwise left rotation of image%s by %g%s",
                                gmic_selection.data(),value,ssep,unsigned int,
                                rol,
                                "Compute bitwise left rotation of image%s by image [%d]",
                                gmic_selection.data(),ind[0],
                                rol,
                                "Compute bitwise left rotation of image%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute sequential bitwise left rotation of image%s");

        // Bitwise right rotation.
        gmic_arithmetic_command("ror",
                                ror,
                                "Compute bitwise right rotation of image%s by %g%s",
                                gmic_selection.data(),value,ssep,unsigned int,
                                ror,
                                "Compute bitwise right rotation of image%s by image [%d]",
                                gmic_selection.data(),ind[0],
                                ror,
                                "Compute bitwise left rotation of image%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute sequential bitwise left rotation of image%s");

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 's...'
        //-----------------------------
      gmic_commands_s :

        // Set status.
        if (!is_get && !std::strcmp("status",item)) {
          gmic_substitute_args(false);
          print(images,0,"Set status to '%s'.",gmic_argument_text_printed());
          CImg<char>::string(argument).move_to(status);
          ++position;
          continue;
        }

        // Skip argument.
        if (is_command_skip) {
          gmic_substitute_args(false);
          if (is_very_verbose)
            print(images,0,"Skip argument '%s'.",
                  gmic_argument_text_printed());
          ++position;
          continue;
        }

        // Set pixel value.
        if (!std::strcmp("set",command)) {
          gmic_substitute_args(false);
          double x = 0, y = 0, z = 0, c = 0;
          value = 0;
          sepx = sepy = sepz = sepc = *argx = *argy = *argz = *argc = 0;
          if ((cimg_sscanf(argument,"%lf%c",
                           &value,&end)==1 ||
               cimg_sscanf(argument,"%lf,%255[0-9.eE%+-]%c",
                           &value,gmic_use_argx,&end)==2 ||
               cimg_sscanf(argument,"%lf,%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           &value,argx,gmic_use_argy,&end)==3 ||
               cimg_sscanf(argument,"%lf,%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           &value,argx,argy,gmic_use_argz,&end)==4 ||
               cimg_sscanf(argument,"%lf,%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-]%c",
                           &value,argx,argy,argz,gmic_use_argc,&end)==5) &&
              (!*argx ||
               (cimg_sscanf(argx,"%lf%c%c",&x,&sepx,&end)==2 && sepx=='%') ||
               cimg_sscanf(argx,"%lf%c",&x,&end)==1) &&
              (!*argy ||
               (cimg_sscanf(argy,"%lf%c%c",&y,&sepy,&end)==2 && sepy=='%') ||
               cimg_sscanf(argy,"%lf%c",&y,&end)==1) &&
              (!*argz ||
               (cimg_sscanf(argz,"%lf%c%c",&z,&sepz,&end)==2 && sepz=='%') ||
               cimg_sscanf(argz,"%lf%c",&z,&end)==1) &&
              (!*argc ||
               (cimg_sscanf(argc,"%lf%c%c",&c,&sepc,&end)==2 && sepc=='%') ||
               cimg_sscanf(argc,"%lf%c",&c,&end)==1)) {
            print(images,0,"Set value %g in image%s, at coordinates (%g%s,%g%s,%g%s,%g%s).",
                  value,
                  gmic_selection.data(),
                  x,sepx=='%'?"%":"",
                  y,sepy=='%'?"%":"",
                  z,sepz=='%'?"%":"",
                  c,sepc=='%'?"%":"");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const int
                nx = (int)cimg::round(sepx=='%'?x*(img.width() - 1)/100:x),
                ny = (int)cimg::round(sepy=='%'?y*(img.height() - 1)/100:y),
                nz = (int)cimg::round(sepz=='%'?z*(img.depth() - 1)/100:z),
                nc = (int)cimg::round(sepc=='%'?c*(img.spectrum() - 1)/100:c);
              gmic_apply(gmic_set(value,nx,ny,nz,nc),true);
            }
          } else arg_error("set");
          is_change = true;
          ++position;
          continue;
        }

        // Split.
        if (!std::strcmp("split",command)) {
          gmic_substitute_args(false);
          double nb = -1;
          char pm = 0;
          *argx = 0;
          if (cimg_sscanf(argument,"%255[xyzc],%lf%c",gmic_use_argx,&nb,&end)==2 ||
              (nb = -1,cimg_sscanf(argument,"%255[xyzc]%c",argx,&end))==1) {

            // Split along axes.
            nb = cimg::round(nb);
            if (nb>0)
              print(images,0,"Split image%s along the '%s'-ax%cs, into %g parts.",
                    gmic_selection.data(),
                    argx,
                    std::strlen(argx)>1?'e':'i',
                    nb);
            else if (nb<0) {
              if (nb==-1)
                print(images,0,"Split image%s along the '%s'-ax%cs.",
                      gmic_selection.data(),
                      argx,
                      std::strlen(argx)>1?'e':'i');
              else
                print(images,0,"Split image%s along the '%s'-ax%cs, into blocs of %g pixels.",
                      gmic_selection.data(),
                      argx,
                      std::strlen(argx)>1?'e':'i',
                      -nb);
            } else
              print(images,0,"Split image%s along the '%s'-ax%cs, according to constant values.",
                    gmic_selection.data(),
                    argx,
                    std::strlen(argx)>1?'e':'i');

            int off = 0;
            cimg_forY(selection,l) {
              uind = selection[l] + off;
              const CImg<T>& img = gmic_check(images[uind]);
              if (!img) {
                if (!is_get) {
                  images.remove(uind);
                  images_names.remove(uind);
                  --off;
                }
              } else {
                g_list.assign(img,true);
                for (const char *p_axis = argx; *p_axis; ++p_axis) {
                  const unsigned int N = g_list.size();
                  for (unsigned int q = 0; q<N; ++q) {
                    g_list[0].get_split(*p_axis,(int)nb).move_to(g_list,~0U);
                    g_list.remove(0);
                  }
                }
                if (is_get) {
                  if (g_list) {
                    pattern = images_names.size();
                    images_names.insert(g_list.size());
                    images_names[uind].get_copymark().move_to(images_names[pattern]);
                    for (unsigned int i = 1; i<g_list.size(); ++i)
                      images_names[pattern + i - 1].get_copymark().move_to(images_names[pattern + i]);
                    g_list.move_to(images,~0U);
                  }
                } else {
                  if (g_list) {
                    images.insert(g_list.size() - 1,uind + 1);
                    images_names.insert(g_list.size() - 1,uind + 1);
                    g_list[0].move_to(images[uind]);
                    for (unsigned int i = 1; i<g_list.size(); ++i) {
                      g_list[i].move_to(images[uind + i]);
                      images_names[uind + i - 1].get_copymark().move_to(images_names[uind + i]);
                    }
                    off+=(int)g_list.size() - 1;
                  } else {
                    images.remove(uind);
                    images_names.remove(uind);
                    --off;
                  }
                }
              }
            }
            g_list.assign();
            is_change = true;
            ++position;
            continue;

          } else if (cimg_sscanf(argument,"%c%c",&pm,&end)==2 && (pm=='+' || pm=='-') && end==',') {

            // Split according to values sequence (opt. with axes too).
            const char *s_values = argument + 2;
            *argx = 0;
            if (cimg_sscanf(s_values,"%255[xyzc]%c",gmic_use_argx,&end)==2 && end==',') s_values+=std::strlen(argx) + 1;
            unsigned int nb_values = 1;
            for (const char *s = s_values; *s; ++s) if (*s==',') ++nb_values;
            CImg<T> values;
            try { values.assign(nb_values,1,1,1).fill_from_values(s_values,true); }
            catch (CImgException&) { values.assign(); }
            if (values) {
              if (*argx)
                print(images,0,"Split image%s in '%s' mode along '%s'-ax%cs, according to value sequence '%s'.",
                      gmic_selection.data(),
                      pm=='-'?"discard":"keep",
                      argx,
                      std::strlen(argx)>1?'e':'i',
                      gmic_argument_text_printed() + (s_values - argument));
              else
                print(images,0,"Split image%s in '%s' mode, according to value sequence '%s'.",
                      gmic_selection.data(),
                      pm=='-'?"discard":"keep",
                      gmic_argument_text_printed() + 2);

              int off = 0;
              cimg_forY(selection,l) {
                uind = selection[l] + off;
                const CImg<T>& img = gmic_check(images[uind]);
                if (!img) {
                  if (!is_get) {
                    images.remove(uind);
                    images_names.remove(uind);
                    --off;
                  }
                } else {
                  if (*argx) { // Along axes
                    g_list.assign(img,true);
                    for (const char *p_axis = argx; *p_axis; ++p_axis) {
                      const unsigned int N = g_list.size();
                      for (unsigned int q = 0; q<N; ++q) {
                        g_list[0].get_split(values,*p_axis,pm=='+').move_to(g_list,~0U);
                        g_list.remove(0);
                      }
                    }
                  } else // Without axes
                    img.get_split(values,0,pm=='+').move_to(g_list);

                  if (is_get) {
                    if (g_list) {
                      pattern = images_names.size();
                      images_names.insert(g_list.size());
                      images_names[uind].get_copymark().move_to(images_names[pattern]);
                      for (unsigned int i = 1; i<g_list.size(); ++i)
                        images_names[pattern + i - 1].get_copymark().move_to(images_names[pattern + i]);
                      g_list.move_to(images,~0U);
                    }
                  } else {
                    if (g_list) {
                      images.insert(g_list.size() - 1,uind + 1);
                      images_names.insert(g_list.size() - 1,uind + 1);
                      g_list[0].move_to(images[uind]);
                      for (unsigned int i = 1; i<g_list.size(); ++i) {
                        g_list[i].move_to(images[uind + i]);
                        images_names[uind + i - 1].get_copymark().move_to(images_names[uind + i]);
                      }
                      off+=(int)g_list.size() - 1;
                    } else {
                      images.remove(uind);
                      images_names.remove(uind);
                      --off;
                    }
                  }
                }
              }
              g_list.assign();
              is_change = true;
              ++position;
              continue;
            }
          }

          // Split as constant one-column vectors.
          print(images,0,"Split image%s as constant one-column vectors.",
                gmic_selection.data());

          int off = 0;
          cimg_forY(selection,l) {
            uind = selection[l] + off;
            CImg<T>& img = gmic_check(images[uind]);
            if (!img) {
              if (!is_get) {
                images.remove(uind);
                images_names.remove(uind);
                --off;
              }
            } else {
              CImg<T>(img.data(),1,(unsigned int)img.size(),1,1,true).get_split('y',0).move_to(g_list);
              if (is_get) {
                if (g_list) {
                  pattern = images_names.size();
                  images_names.insert(g_list.size());
                  images_names[uind].get_copymark().move_to(images_names[pattern]);
                  for (unsigned int i = 1; i<g_list.size(); ++i)
                    images_names[pattern + i - 1].get_copymark().move_to(images_names[pattern + i]);
                  g_list.move_to(images,~0U);
                }
              } else {
                if (g_list) {
                  images.insert(g_list.size() - 1,uind + 1);
                  images_names.insert(g_list.size() - 1,uind + 1);
                  g_list[0].move_to(images[uind]);
                  for (unsigned int i = 1; i<g_list.size(); ++i) {
                    g_list[i].move_to(images[uind + i]);
                    images_names[uind + i - 1].get_copymark().move_to(images_names[uind + i]);
                  }
                  off+=(int)g_list.size() - 1;
                } else {
                  images.remove(uind);
                  images_names.remove(uind);
                  --off;
                }
              }
            }
          }
          g_list.assign();
          is_change = true;
          continue;
        }

        // Shared input.
        if (!std::strcmp("shared",command)) {
          gmic_substitute_args(false);
          CImg<char> st0(256), st1(256), st2(256), st3(256), st4(256);
          char sep2 = 0, sep3 = 0, sep4 = 0;
          double a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0;
          sep0 = sep1 = *st0 = *st1 = *st2 = *st3 = *st4 = 0;
          pattern = images.size();
          images.insert(selection.height());
          images_names.insert(selection.height());
          if (cimg_sscanf(argument,
                          "%255[0-9.eE%+],%255[0-9.eE%+],%255[0-9.eE%+],%255[0-9.eE%+],"
                          "%255[0-9.eE%+]%c",
                          st0.data(),st1.data(),st2.data(),st3.data(),st4.data(),&end)==5 &&
              (cimg_sscanf(st0,"%lf%c",&a0,&end)==1 ||
               (cimg_sscanf(st0,"%lf%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
              (cimg_sscanf(st1,"%lf%c",&a1,&end)==1 ||
               (cimg_sscanf(st1,"%lf%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
              (cimg_sscanf(st2,"%lf%c",&a2,&end)==1 ||
               (cimg_sscanf(st2,"%lf%c%c",&a2,&sep2,&end)==2 && sep2=='%')) &&
              (cimg_sscanf(st3,"%lf%c",&a3,&end)==1 ||
               (cimg_sscanf(st3,"%lf%c%c",&a3,&sep3,&end)==2 && sep3=='%')) &&
              (cimg_sscanf(st4,"%lf%c",&a4,&end)==1 ||
               (cimg_sscanf(st4,"%lf%c%c",&a4,&sep4,&end)==2 && sep4=='%'))) {
            print(images,0,
                  "Insert shared buffer%s from points (%g%s->%g%s,%g%s,%g%s,%g%s) of image%s.",
                  selection.height()>1?"s":"",
                  a0,sep0=='%'?"%":"",
                  a1,sep1=='%'?"%":"",
                  a2,sep2=='%'?"%":"",
                  a3,sep3=='%'?"%":"",
                  a4,sep4=='%'?"%":"",
                  gmic_selection.data());

            cimg_forY(selection,l) {
              CImg<T>& img = images[selection[l]];
              nvalue0 = cimg::round(sep0=='%'?a0*(img.width() - 1)/100:a0);
              nvalue1 = cimg::round(sep1=='%'?a1*(img.width() - 1)/100:a1);
              const unsigned int
                y =  (unsigned int)cimg::round(sep2=='%'?a2*(img.height() - 1)/100:a2),
                z =  (unsigned int)cimg::round(sep3=='%'?a3*(img.depth() - 1)/100:a3),
                c =  (unsigned int)cimg::round(sep4=='%'?a4*(img.spectrum() - 1)/100:a4);
              images[pattern + l].assign(img.get_shared_points((unsigned int)nvalue0,(unsigned int)nvalue1,y,z,c),true);
              images_names[selection[l]].get_copymark().move_to(images_names[pattern + l]);
            }
            ++position;
          } else if (cimg_sscanf(argument,
                                 "%255[0-9.eE%+],%255[0-9.eE%+],%255[0-9.eE%+],"
                                 "%255[0-9.eE%+],%c",
                                 st0.data(),st1.data(),st2.data(),st3.data(),&end)==4 &&
                     (cimg_sscanf(st0,"%lf%c",&a0,&end)==1 ||
                      (cimg_sscanf(st0,"%lf%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (cimg_sscanf(st1,"%lf%c",&a1,&end)==1 ||
                      (cimg_sscanf(st1,"%lf%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
                     (cimg_sscanf(st2,"%lf%c",&a2,&end)==1 ||
                      (cimg_sscanf(st2,"%lf%c%c",&a2,&sep2,&end)==2 && sep2=='%')) &&
                     (cimg_sscanf(st3,"%lf%c",&a3,&end)==1 ||
                      (cimg_sscanf(st3,"%lf%c%c",&a3,&sep3,&end)==2 && sep3=='%'))) {
            print(images,0,"Insert shared buffer%s from lines (%g%s->%g%s,%g%s,%g%s) of image%s.",
                  selection.height()>1?"s":"",
                  a0,sep0=='%'?"%":"",
                  a1,sep1=='%'?"%":"",
                  a2,sep2=='%'?"%":"",
                  a3,sep3=='%'?"%":"",
                  gmic_selection.data());

            cimg_forY(selection,l) {
              CImg<T>& img = images[selection[l]];
              nvalue0 = cimg::round(sep0=='%'?a0*(img.height() - 1)/100:a0);
              nvalue1 = cimg::round(sep1=='%'?a1*(img.height() - 1)/100:a1);
              const unsigned int
                z =  (unsigned int)cimg::round(sep2=='%'?a2*(img.depth() - 1)/100:a2),
                c =  (unsigned int)cimg::round(sep3=='%'?a3*(img.spectrum() - 1)/100:a3);
              images[pattern + l].assign(img.get_shared_rows((unsigned int)nvalue0,(unsigned int)nvalue1,z,c),true);
              images_names[selection[l]].get_copymark().move_to(images_names[pattern + l]);
            }
            ++position;
          } else if (cimg_sscanf(argument,"%255[0-9.eE%+],%255[0-9.eE%+],%255[0-9.eE%+]%c",
                                 st0.data(),st1.data(),st2.data(),&end)==3 &&
                     (cimg_sscanf(st0,"%lf%c",&a0,&end)==1 ||
                      (cimg_sscanf(st0,"%lf%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (cimg_sscanf(st1,"%lf%c",&a1,&end)==1 ||
                      (cimg_sscanf(st1,"%lf%c%c",&a1,&sep1,&end)==2 && sep1=='%')) &&
                     (cimg_sscanf(st2,"%lf%c",&a2,&end)==1 ||
                      (cimg_sscanf(st2,"%lf%c%c",&a2,&sep2,&end)==2 && sep2=='%'))) {
            print(images,0,"Insert shared buffer%s from planes (%g%s->%g%s,%g%s) of image%s.",
                  selection.height()>1?"s":"",
                  a0,sep0=='%'?"%":"",
                  a1,sep1=='%'?"%":"",
                  a2,sep2=='%'?"%":"",
                  gmic_selection.data());
            cimg_forY(selection,l) {
              CImg<T>& img = images[selection[l]];
              nvalue0 = cimg::round(sep0=='%'?a0*(img.depth() - 1)/100:a0);
              nvalue1 = cimg::round(sep1=='%'?a1*(img.depth() - 1)/100:a1);
              const unsigned int
                c =  (unsigned int)cimg::round(sep2=='%'?a2*(img.spectrum() - 1)/100:a2);
              images[pattern + l].assign(img.get_shared_slices((unsigned int)nvalue0,(unsigned int)nvalue1,c),true);
              images_names[selection[l]].get_copymark().move_to(images_names[pattern + l]);
            }
            ++position;
          } else if (cimg_sscanf(argument,"%255[0-9.eE%+],%255[0-9.eE%+]%c",
                                 st0.data(),st1.data(),&end)==2 &&
                     (cimg_sscanf(st0,"%lf%c",&a0,&end)==1 ||
                      (cimg_sscanf(st0,"%lf%c%c",&a0,&sep0,&end)==2 && sep0=='%')) &&
                     (cimg_sscanf(st1,"%lf%c",&a1,&end)==1 ||
                      (cimg_sscanf(st1,"%lf%c%c",&a1,&sep1,&end)==2 && sep1=='%'))) {
            print(images,0,"Insert shared buffer%s from channels (%g%s->%g%s) of image%s.",
                  selection.height()>1?"s":"",
                  a0,sep0=='%'?"%":"",
                  a1,sep1=='%'?"%":"",
                  gmic_selection.data());
            cimg_forY(selection,l) {
              CImg<T>& img = images[selection[l]];
              nvalue0 = cimg::round(sep0=='%'?a0*(img.spectrum() - 1)/100:a0);
              nvalue1 = cimg::round(sep1=='%'?a1*(img.spectrum() - 1)/100:a1);
              images[pattern + l].assign(img.get_shared_channels((unsigned int)nvalue0,(unsigned int)nvalue1),true);
              images_names[selection[l]].get_copymark().move_to(images_names[pattern + l]);
            }
            ++position;
          } else if (cimg_sscanf(argument,"%255[0-9.eE%+]%c",
                                 st0.data(),&end)==1 &&
                     (cimg_sscanf(st0,"%lf%c",&a0,&end)==1 ||
                      (cimg_sscanf(st0,"%lf%c%c",&a0,&sep0,&end)==2 && sep0=='%'))) {
            print(images,0,"Insert shared buffer%s from channel %g%s of image%s.",
                  selection.height()>1?"s":"",
                  a0,sep0=='%'?"%":"",
                  gmic_selection.data());
            cimg_forY(selection,l) {
              CImg<T>& img = images[selection[l]];
              nvalue0 = cimg::round(sep0=='%'?a0*(img.spectrum() - 1)/100:a0);
              images[pattern + l].assign(img.get_shared_channel((unsigned int)nvalue0),true);
              images_names[selection[l]].get_copymark().move_to(images_names[pattern + l]);
            }
            ++position;
          } else {
            print(images,0,"Insert shared buffer%s from image%s.",
                  selection.height()>1?"s":"",
                  gmic_selection.data());
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              images[pattern + l].assign(img,true);
              images_names[selection[l]].get_copymark().move_to(images_names[pattern + l]);
            }
          }
          is_change = true;
          continue;
        }

        // Store.
        if (!std::strcmp("store",command)) {
          gmic_substitute_args(false);
          unsigned int is_compressed = 0U;
          if ((cimg_sscanf(argument,"%u,%4095[,a-zA-Z0-9_]%c",&is_compressed,&(*gmic_use_formula=0),&end)==2 ||
               cimg_sscanf(argument,"%4095[,a-zA-Z0-9_]%c",&(*formula=0),&end)==1) &&
              is_compressed<=1 &&
              (*formula<'0' || *formula>'9') && *formula!=',') {
            char *current = formula, *next = std::strchr(current,','), saved = 0;
            pattern = 1U;
            if (next) // Count number of specified variable names
              for (const char *ptr = next; ptr; ptr = std::strchr(ptr,',')) { ++ptr; ++pattern; }
            print(images,0,
                  "Store image%s as %svariable%s '%s'",
                  gmic_selection.data(),
                  is_compressed?"compressed":"",
                  next?"s":"",
                  gmic_argument_text_printed() + (*argument=='0' || *argument=='1'?2:0));

            if (pattern!=1 && (int)pattern!=selection.height())
              error(true,images,0,0,
                    "Command 'store': Specified arguments '%s' do not match numbers of selected images.",
                    gmic_argument_text() + (*argument=='0' || *argument=='1'?2:0));

            g_list.assign(selection.height());
            g_list_c.assign(g_list.size());
            cimg_forY(selection,l) {
              uind = selection[l];
              if (is_get) {
                g_list[l] = images[uind].get_shared();
                CImg<char>::string(images_names[uind]).move_to(g_list_c[l]);
              } else {
                images[uind].move_to(g_list[l]);
                CImg<char>::string(images_names[uind]).move_to(g_list_c[l]);
                images_names[uind].assign();
              }
            }

            if (pattern==1) { // Assignment to a single variable
              (g_list_c>'x').move_to(name);
              name.resize(name.width() + 4,1,1,1,0,0,1);
              name[0] = 'G'; name[1] = 'M'; name[2] = 'Z'; name[3] = 0;
              name.unroll('y').move_to(g_list);
              g_list.get_serialize((bool)is_compressed,(unsigned int)(9 + std::strlen(formula))).move_to(name);
              cimg_snprintf(name,name._height,"%c*store/%s",gmic_store,_formula.data());
              set_variable(formula,name,variables_sizes);
            } else for (unsigned int n = 0; n<pattern; ++n) { // Assignment to multiple variables
              if (!*current)
                error(true,images,0,0,
                      "Command 'store': Empty variable name specified in arguments '%s'.",
                      gmic_argument_text() + (*argument=='0' || *argument=='1'?2:0));
              saved = next?*next:0;
              if (next) *next = 0;

              CImgList<T> tmp(2);
              g_list_c[n].move_to(name);
              name.resize(name.width() + 4,1,1,1,0,0,1);
              name[0] = 'G'; name[1] = 'M'; name[2] = 'Z'; name[3] = 0;
              name.unroll('y').move_to(tmp[1]);
              g_list[n].move_to(tmp[0]);
              tmp.get_serialize((bool)is_compressed,(unsigned int)(9 + std::strlen(current))).move_to(name);
              cimg_snprintf(name,name._height,"%c*store/%s",gmic_store,current);
              set_variable(current,name,variables_sizes);

              if (saved) { // Other variables names follow
                *next = saved;
                current = next + 1;
                next = std::strchr(next + 1,',');
              }
            }
            if (!is_get) remove_images(images,images_names,selection,0,selection.height() - 1);
            name.assign();
            g_list.assign();
            g_list_c.assign();
          } else arg_error("store");
          ++position;
          continue;
        }

        // Shift.
        if (!std::strcmp("shift",command)) {
          gmic_substitute_args(false);
          double dx = 0, dy = 0, dz = 0, dc = 0;
          sepx = sepy = sepz = sepc = *argx = *argy = *argz = *argc = 0;
          interpolation = boundary = 0;
          if ((cimg_sscanf(argument,"%255[0-9.eE%+-]%c",
                           gmic_use_argx,&end)==1 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           argx,gmic_use_argy,&end)==2 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           argx,argy,gmic_use_argz,&end)==3 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-]%c",
                           argx,argy,argz,gmic_use_argc,&end)==4 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-],%u%c",
                           argx,argy,argz,argc,&boundary,&end)==5 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-],%u,%u%c",
                           argx,argy,argz,argc,&boundary,&interpolation,&end)==6) &&
              (cimg_sscanf(argx,"%lf%c",&dx,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&dx,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy ||
               cimg_sscanf(argy,"%lf%c",&dy,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&dy,&sepy,&end)==2 && sepy=='%')) &&
              (!*argz ||
               cimg_sscanf(argz,"%lf%c",&dz,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&dz,&sepz,&end)==2 && sepz=='%')) &&
              (!*argc ||
               cimg_sscanf(argc,"%lf%c",&dc,&end)==1 ||
               (cimg_sscanf(argc,"%lf%c%c",&dc,&sepc,&end)==2 && sepc=='%')) &&
              boundary<=3 && interpolation<=1) {
            print(images,0,
                  "Shift image%s by displacement vector (%g%s,%g%s,%g%s,%g%s), "
                  "%s boundary conditions and %s interpolation.",
                  gmic_selection.data(),
                  dx,sepx=='%'?"%":"",
                  dy,sepy=='%'?"%":"",
                  dz,sepz=='%'?"%":"",
                  dc,sepc=='%'?"%":"",
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror",
                  interpolation?"linear":"nearest-neighbor");
            cimg_forY(selection,l) {
              CImg<T> &img = images[selection[l]];
              const float
                ndx = (float)cimg::round(sepx=='%'?dx*img.width()/100:dx),
                ndy = (float)cimg::round(sepy=='%'?dy*img.height()/100:dy),
                ndz = (float)cimg::round(sepz=='%'?dz*img.depth()/100:dz),
                ndc = (float)cimg::round(sepc=='%'?dc*img.spectrum()/100:dc);
              gmic_apply(gmic_shift(ndx,ndy,ndz,ndc,boundary,(bool)interpolation),false);
            }
          } else arg_error("shift");
          is_change = true;
          ++position;
          continue;
        }

        // Sub.
        gmic_arithmetic_command("sub",
                                operator-=,
                                "Subtract %g%s to image%s",
                                value,ssep,gmic_selection.data(),Tfloat,
                                operator-=,
                                "Subtract image [%d] to image%s",
                                ind[0],gmic_selection.data(),
                                operator_minuseq,
                                "Subtract expression %s to image%s",
                                gmic_argument_text_printed(),gmic_selection.data(),
                                "Subtract image%s");
        // Square root.
        gmic_simple_command("sqrt",sqrt,"Compute pointwise square root of image%s.");

        // Square.
        gmic_simple_command("sqr",sqr,"Compute pointwise square function of image%s.");

        // Sign.
        gmic_simple_command("sign",sign,"Compute pointwise sign of image%s.");

        // Sine.
        gmic_simple_command("sin",sin,"Compute pointwise sine of image%s.");

        // Sort.
        if (!std::strcmp("sort",command)) {
          gmic_substitute_args(false);
          char order = '+';
          axis = 0;
          if ((cimg_sscanf(argument,"%c%c",&order,&end)==1 ||
               (cimg_sscanf(argument,"%c,%c%c",&order,&axis,&end)==2 &&
                is_xyzc(axis))) && (order=='+' || order=='-')) ++position;
          else { order = '+'; axis = 0; }
          if (axis) print(images,0,"Sort values of image%s in %s order, according to axis '%c'.",
                          gmic_selection.data(),order=='+'?"ascending":"descending",axis);
          else print(images,0,"Sort values of image%s in %s order.",
                     gmic_selection.data(),order=='+'?"ascending":"descending");
          cimg_forY(selection,l) gmic_apply(sort(order=='+',axis),false);
          is_change = true;
          continue;
        }

        // Solve.
        if (!std::strcmp("solve",command)) {
          gmic_substitute_args(true);
          sep = *indices = 0;
          pattern = 0;
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep,&end)==2 && sep==']') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u%c",gmic_use_indices,&pattern,&end)==2) &&
              (ind=selection2cimg(indices,images.size(),images_names,"solve")).height()==1 &&
              pattern<=1) {
            print(images,0,"Solve linear system AX = B, with B-vector%s, A-matrix [%d] and %s solver.",
                  gmic_selection.data(),*ind,pattern?"LU":"SVD");
            const CImg<double> A = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply_double(solve(A,(bool)pattern));
          } else arg_error("solve");
          is_change = true;
          ++position;
          continue;
        }

        // Shift 3D object, with opposite displacement.
        if (!std::strcmp("sub3d",command)) {
          gmic_substitute_args(false);
          float tx = 0, ty = 0, tz = 0;
          if (cimg_sscanf(argument,"%f%c",
                          &tx,&end)==1 ||
              cimg_sscanf(argument,"%f,%f%c",
                          &tx,&ty,&end)==2 ||
              cimg_sscanf(argument,"%f,%f,%f%c",
                          &tx,&ty,&tz,&end)==3) {
            print(images,0,"Shift 3D object%s with displacement -(%g,%g,%g).",
                  gmic_selection.data(),
                  tx,ty,tz);
            cimg_forY(selection,l) {
              uind = selection[l];
              CImg<T>& img = images[uind];
              try { gmic_apply(shift_CImg3d(-tx,-ty,-tz),true); }
              catch (CImgException&) {
                if (!img.is_CImg3d(true,&(*gmic_use_message=0)))
                  error(true,images,0,0,
                        "Command 'sub3d': Invalid 3D object [%d], in image%s (%s).",
                        uind,gmic_selection_err.data(),message);
                else throw;
              }
            }
          } else arg_error("sub3d");
          is_change = true;
          ++position;
          continue;
        }

        // Set random generator seed.
        if (!is_get && !std::strcmp("srand",item)) {
          gmic_substitute_args(false);
          value = 0;
          if (cimg_sscanf(argument,"%lf%c",
                          &value,&end)==1) {
            value = cimg::round(value);
            print(images,0,"Set random generator seed to %u.",
                  (unsigned int)value);
            cimg::srand((unsigned int)value);
            ++position;
          } else {
            print(images,0,"Set random generator seed to random.");
            cimg::srand();
          }
          continue;
        }

        // Anisotropic PDE-based smoothing.
        if (!std::strcmp("smooth",command)) {
          gmic_substitute_args(true);
          float sharpness = 0.7f, anisotropy = 0.3f, dl =0.8f, da = 30.f, gauss_prec = 2.f;
          unsigned int is_fast_approximation = 1;
          *argx = *argy = *argz = sep = sep0 = sep1 = 0;
          interpolation = 0;
          value0 = 0.6; value1 = 1.1;
          if ((cimg_sscanf(argument,"%255[0-9.eE%+-]%c",
                           gmic_use_argz,&end)==1 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%f%c",
                           argz,&sharpness,&end)==2 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%f,%f%c",
                           argz,&sharpness,&anisotropy,&end)==3 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%f,%f,%255[0-9.eE%+-]%c",
                           argz,&sharpness,&anisotropy,gmic_use_argx,&end)==4 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%f,%f,%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           argz,&sharpness,&anisotropy,argx,gmic_use_argy,&end)==5 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%f,%f,%255[0-9.eE%+-],%255[0-9.eE%+-],%f%c",
                           argz,&sharpness,&anisotropy,argx,argy,&dl,&end)==6 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%f,%f,%255[0-9.eE%+-],%255[0-9.eE%+-],%f,%f%c",
                           argz,&sharpness,&anisotropy,argx,argy,&dl,&da,&end)==7 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%f,%f,%255[0-9.eE%+-],%255[0-9.eE%+-],%f,%f,%f%c",
                           argz,&sharpness,&anisotropy,argx,argy,&dl,&da,&gauss_prec,&end)==8 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%f,%f,%255[0-9.eE%+-],%255[0-9.eE%+-],%f,%f,%f,%u%c",
                           argz,&sharpness,&anisotropy,argx,argy,&dl,&da,&gauss_prec,&interpolation,&end)==9 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%f,%f,%255[0-9.eE%+-],%255[0-9.eE%+-],%f,%f,%f,%u,%u,%c",
                           argz,&sharpness,&anisotropy,argx,argy,&dl,&da,&gauss_prec,&interpolation,
                           &is_fast_approximation,&end)==10) &&
              (cimg_sscanf(argz,"%lf%c",&value,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&value,&sep,&end)==2 && sep=='%')) &&
              (!*argx ||
               cimg_sscanf(argx,"%lf%c",&value0,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&value0,&sep0,&end)==2 && sep0=='%')) &&
              (!*argy ||
               cimg_sscanf(argy,"%lf%c",&value1,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&value1,&sep1,&end)==2 && sep1=='%')) &&
              value>=0 && value0>=0 && value1>=0 && sharpness>=0 && anisotropy>=0 && anisotropy<=1 && dl>0 &&
              (da>0 || (da==0 && sep!='%')) && gauss_prec>0 && interpolation<=2 && is_fast_approximation<=1) {
            if (da>0)
              print(images,0,"Smooth image%s anisotropically, with amplitude %g%s, sharpness %g, "
                    "anisotropy %g, alpha %g%s, sigma %g%s, dl %g, da %g, precision %g, "
                    "%s interpolation and fast approximation %s.",
                    gmic_selection.data(),
                    value,sep=='%'?"%":"",
                    sharpness,anisotropy,value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",dl,da,gauss_prec,
                    interpolation==0?"nearest-neighbor":interpolation==1?"linear":"runge-kutta",
                    is_fast_approximation?"enabled":"disabled");
            else {
              value = cimg::round(value);
              print(images,0,"Smooth image%s anisotropically, with %d iterations, sharpness %g, "
                    "anisotropy %g, alpha %g%s, sigma %g%s and dt %g.",
                    gmic_selection.data(),
                    (int)value,
                    sharpness,anisotropy,value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",dl);
            }
            if (sep=='%') value = -value;
            if (sep0=='%') value0 = -value0;
            if (sep1=='%') value1 = -value1;
            cimg_forY(selection,l)
              gmic_apply(blur_anisotropic((float)value,sharpness,anisotropy,(float)value0,(float)value1,dl,da,
                                          gauss_prec,interpolation,(bool)is_fast_approximation),true);
          } else if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                                   gmic_use_indices,&sep,&end)==2 && sep==']') ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-]%c",
                                  indices,gmic_use_argx,&end)==2 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%f%c",
                                  indices,argx,&dl,&end)==3 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%f,%f%c",
                                  indices,argx,&dl,&da,&end)==4 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%f,%f,%f%c",
                                  indices,argx,&dl,&da,&gauss_prec,&end)==5 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%f,%f,%f,%u%c",
                                  indices,argx,&dl,&da,&gauss_prec,&interpolation,&end)==6 ||
                      cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%255[0-9.eE%+-],%f,%f,%f,%u,%u%c",
                                  indices,argx,&dl,&da,&gauss_prec,&interpolation,
                                  &is_fast_approximation,&end)==7) &&
                     (ind=selection2cimg(indices,images.size(),images_names,"smooth")).height()==1 &&
                     (!*argx ||
                      cimg_sscanf(argx,"%lf%c",&value0,&end)==1 ||
                      (cimg_sscanf(argx,"%lf%c%c",&value0,&sep0,&end)==2 && sep0=='%')) &&
                     value0>=0 && dl>0 && (da>0 || (da==0 && sep0!='%')) && gauss_prec>0 && interpolation<=2 &&
                     is_fast_approximation<=1) {
            const CImg<T> tensors = gmic_image_arg(*ind);
            if (da>0)
              print(images,0,
                    "Smooth image%s anisotropically, with tensor field [%u], amplitude %g%s, "
                    "dl %g, da %g, precision %g, %s interpolation and fast approximation %s.",
                    gmic_selection.data(),
                    *ind,
                    value0,sep0=='%'?"%":"",
                    dl,da,gauss_prec,interpolation==0?"nearest-neighbor":interpolation==1?"linear":"runge-kutta",
                    is_fast_approximation?"enabled":"disabled");
            else {
              value0 = cimg::round(value0);
              print(images,0,
                    "Smooth image%s anisotropically, with tensor field [%u], %d iterations "
                    "and dt %g.",
                    gmic_selection.data(),
                    *ind,(int)value0,dl);
            }
            cimg_forY(selection,l)
              gmic_apply(blur_anisotropic(tensors,(float)value0,dl,da,gauss_prec,interpolation,
                                          is_fast_approximation),true);
          } else arg_error("smooth");
          is_change = true;
          ++position;
          continue;
        }

        // Screenshot.
        if (!std::strcmp("screen",item)) {
          gmic_substitute_args(false);
          gmic_use_title;
          sepx = sepy = sepz = sepc = *argx = *argy = *argz = *argc = 0;
          if (cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                          gmic_use_argx,gmic_use_argy,gmic_use_argz,gmic_use_argc,&end)==4 &&
              (cimg_sscanf(argx,"%lf%c",&value0,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&value0,&sepx,&end)==2 && sepx=='%')) &&
              (cimg_sscanf(argy,"%lf%c",&value1,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&value1,&sepy,&end)==2 && sepy=='%')) &&
              (cimg_sscanf(argz,"%lf%c",&nvalue0,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&nvalue0,&sepz,&end)==2 && sepz=='%')) &&
              (cimg_sscanf(argc,"%lf%c",&nvalue1,&end)==1 ||
               (cimg_sscanf(argc,"%lf%c%c",&nvalue1,&sepc,&end)==2 && sepc=='%'))) {
            print(images,0,"Take screenshot, with coordinates (%s,%s) - (%s,%s).",
                  argx,argy,argz,argc);
            if (sepx=='%') value0 = value0*CImgDisplay::screen_width()/100;
            if (sepy=='%') value1 = value1*CImgDisplay::screen_height()/100;
            if (sepz=='%') nvalue0 = nvalue0*CImgDisplay::screen_width()/100;
            if (sepc=='%') nvalue1 = nvalue1*CImgDisplay::screen_height()/100;
            value0 = cimg::round(value0);
            value1 = cimg::round(value1);
            nvalue0 = cimg::round(nvalue0);
            nvalue1 = cimg::round(nvalue1);
            cimg_snprintf(title,_title.width(),"[Screenshot (%g,%g)-(%g,%g)]",
                          value0,value1,nvalue0,nvalue1);
            CImgDisplay::screenshot((int)value0,(int)value1,(int)nvalue0,(int)nvalue1,images.insert(1).back());
            ++position;
          } else {
            print(images,0,"Take screenshot.");
            cimg_snprintf(title,_title.width(),"[Screenshot]");
            CImgDisplay::screenshot(images.insert(1).back());
          }
          CImg<char>::string(title).move_to(images_names);
          is_change = true;
          continue;
        }

        // SVD.
        if (!std::strcmp("svd",command)) {
          print(images,0,"Compute SVD decomposition%s of matri%s%s.",
                selection.height()>1?"s":"",selection.height()>1?"ce":"x",gmic_selection.data());
          CImg<float> U, S, V;
          unsigned int off = 0;
          cimg_forY(selection,l) {
            uind = selection[l] + off;
            const CImg<T>& img = gmic_check(images[uind]);
            img.SVD(U,S,V,true,100);
            if (is_get) {
              U.move_to(images);
              S.move_to(images);
              V.move_to(images);
              images_names[uind].get_copymark().move_to(images_names);
              images_names.back().get_copymark().move_to(images_names);
              images_names.back().get_copymark().move_to(images_names);
            } else {
              images.insert(2,uind + 1);
              U.move_to(images[uind].assign());
              S.move_to(images[uind + 1]);
              V.move_to(images[uind + 2]);
              images_names.insert(2,uind + 1);
              images_names[uind].get_copymark().move_to(images_names[uind + 1]);
              images_names[uind + 1].get_copymark().move_to(images_names[uind + 2]);
              off+=2;
            }
          }
          is_change = true;
          continue;
        }

        // Sine-cardinal.
        gmic_simple_command("sinc",sinc,"Compute pointwise sinc function of image%s.");

        // Hyperbolic sine.
        gmic_simple_command("sinh",sinh,"Compute pointwise hyperpolic sine of image%s.");

        // Extract 3D streamline.
        if (!std::strcmp("streamline3d",command)) {
          gmic_substitute_args(false);
          unsigned int is_backward = 0, is_oriented_only = 0;
          double x = 0, y = 0, z = 0, L = 100, dl = 0.1f;
          sepx = sepy = sepz = *formula = 0;
          interpolation = 2;
          if ((cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           gmic_use_argx,gmic_use_argy,gmic_use_argz,&end)==3 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%lf%c",
                           argx,argy,argz,&L,&end)==4 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%lf,%lf%c",
                           argx,argy,argz,&L,&dl,&end)==5 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%lf,%lf,%u%c",
                           argx,argy,argz,&L,&dl,&interpolation,&end)==6 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%lf,%lf,%u,"
                           "%u%c",
                           argx,argy,argz,&L,&dl,&interpolation,&is_backward,&end)==7 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%lf,%lf,%u,"
                           "%u,%u%c",
                           argx,argy,argz,&L,&dl,&interpolation,&is_backward,
                           &is_oriented_only,&end)==8) &&
              (cimg_sscanf(argx,"%lf%c",&x,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&x,&sepx,&end)==2 && sepx=='%')) &&
              (!*argy ||
               cimg_sscanf(argy,"%lf%c",&y,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&y,&sepy,&end)==2 && sepy=='%')) &&
              (!*argz ||
               cimg_sscanf(argz,"%lf%c",&z,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&z,&sepz,&end)==2 && sepz=='%')) &&
              L>=0 && dl>0 && interpolation<4 && is_backward<=1 && is_oriented_only<=1) {
            print(images,0,"Extract 3D streamline from image%s, starting from (%g%s,%g%s,%g%s).",
                  gmic_selection.data(),
                  x,sepx=='%'?"%":"",
                  y,sepy=='%'?"%":"",
                  z,sepz=='%'?"%":"");
            cimg_forY(selection,l) {
              uind = selection[l];
              CImg<T>& img = gmic_check(images[uind]);
              const float
                nx = (float)cimg::round(sepx=='%'?x*(img.width() - 1)/100:x),
                ny = (float)cimg::round(sepy=='%'?y*(img.height() - 1)/100:y),
                nz = (float)cimg::round(sepz=='%'?z*(img.depth() - 1)/100:z);
              img.get_streamline(nx,ny,nz,(float)L,(float)dl,interpolation,
                                 (bool)is_backward,(bool)is_oriented_only).move_to(vertices);
              if (vertices.width()>1) {
                primitives.assign(vertices.width() - 1,1,2);
                cimglist_for(primitives,q) { primitives(q,0) = (unsigned int)q; primitives(q,1) = q + 1U; }
                g_list_uc.assign(primitives.size(),1,3,1,1,200);
              } else {
                vertices.assign();
                warn(images,0,false,
                     "Command 'streamline3d': Empty streamline starting from "
                     "(%g%s,%g%s,%g%s) in image [%u].",
                     x,sepx=='%'?"%":"",
                     y,sepy=='%'?"%":"",
                     z,sepz=='%'?"%":"",
                     uind);
              }
              vertices.object3dtoCImg3d(primitives,g_list_uc,false);
              gmic_apply(replace(vertices),false);
              primitives.assign();
              g_list_uc.assign();
            }
          } else if ((cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf%c",
                                  gmic_use_formula,&x,&y,&z,&end)==4 ||
                      cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf%c",
                                  formula,&x,&y,&z,&L,&end)==5 ||
                      cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf%c",
                                  formula,&x,&y,&z,&L,&dl,&end)==6 ||
                      cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%u%c",
                                  formula,&x,&y,&z,&L,&dl,&interpolation,&end)==7 ||
                      cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%u,%u%c",
                                  formula,&x,&y,&z,&L,&dl,&interpolation,&is_backward,&end)==8 ||
                      cimg_sscanf(argument,"'%4095[^']',%lf,%lf,%lf,%lf,%lf,%u,%u,%u%c",
                                  formula,&x,&y,&z,&L,&dl,&interpolation,&is_backward,
                                  &is_oriented_only,&end)==9) &&
                     dl>0 && interpolation<4) {
            strreplace_fw(formula);
            print(images,0,"Extract 3D streamline from formula '%s', starting from (%g,%g,%g).",
                  formula,
                  x,y,z);
            CImg<T>::streamline((const char *)formula,(float)x,(float)y,(float)z,(float)L,(float)dl,interpolation,
                                (bool)is_backward,(bool)is_oriented_only).move_to(vertices);
            if (vertices.width()>1) {
              primitives.assign(vertices.width() - 1,1,2);
              cimglist_for(primitives,l) { primitives(l,0) = (unsigned int)l; primitives(l,1) = l + 1U; }
              g_list_uc.assign(primitives.size(),1,3,1,1,200);
            } else {
              vertices.assign();
              warn(images,0,false,
                   "Command 'streamline3d': Empty streamline starting from (%g,%g,%g) "
                   "in expression '%s'.",
                   x,y,z,formula);
            }
            vertices.object3dtoCImg3d(primitives,g_list_uc,false).move_to(images);
            primitives.assign();
            g_list_uc.assign();
            gmic_use_title;
            cimg_snprintf(title,_title.width(),"[3D streamline of '%s' at (%g,%g,%g)]",
                          formula,x,y,z);
            CImg<char>::string(title).move_to(images_names);
          } else arg_error("streamline3d");
          is_change = true;
          ++position;
          continue;
        }

        // Select image feature.
        if (!std::strcmp("select",command)) {
          gmic_substitute_args(false);
          unsigned int feature_type = 0, is_deep_selection = 0;
          *argx = *argy = *argz = sep = sep0 = sep1 = 0;
          value = value0 = value1 = 0;
          exit_on_anykey = 0;
          if ((cimg_sscanf(argument,"%u%c",&feature_type,&end)==1 ||
               (cimg_sscanf(argument,"%u,%255[0-9.eE%+-]%c",
                            &feature_type,gmic_use_argx,&end)==2) ||
               (cimg_sscanf(argument,"%u,%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                            &feature_type,argx,gmic_use_argy,&end)==3) ||
               (cimg_sscanf(argument,"%u,%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                            &feature_type,argx,argy,gmic_use_argz,&end)==4) ||
               (cimg_sscanf(argument,"%u,%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%u%c",
                            &feature_type,argx,argy,argz,&exit_on_anykey,&end)==5) ||
               (cimg_sscanf(argument,"%u,%255[0-9.eE%+-],%255[0-9.eE%+-],%255[0-9.eE%+-],%u,%u%c",
                            &feature_type,argx,argy,argz,&exit_on_anykey,&is_deep_selection,&end)==6)) &&
              (!*argx ||
               cimg_sscanf(argx,"%lf%c",&value,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&value,&sep,&end)==2 && sep=='%')) &&
              (!*argy ||
               cimg_sscanf(argy,"%lf%c",&value0,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&value0,&sep0,&end)==2 && sep0=='%')) &&
              (!*argz ||
               cimg_sscanf(argz,"%lf%c",&value1,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&value1,&sep1,&end)==2 && sep1=='%')) &&
              value>=0 && value0>=0 && value1>=0 && feature_type<=3 && exit_on_anykey<=1 && is_deep_selection<=1) {
            if (!*argx) { value = 50; sep = '%'; }
            if (!*argy) { value0 = 50; sep0 = '%'; }
            if (!*argz) { value1 = 50; sep1 = '%'; }

            if (!is_display_available) {
              print(images,0,
                    "Select %s in image%s in interactive mode, from point (%g%s,%g%s,%g%s) (skipped no display %s).",
                    feature_type==0?"point":feature_type==1?"segment":
                    feature_type==2?"rectangle":"ellipse",gmic_selection.data(),
                    value,sep=='%'?"%":"",value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"",
                    cimg_display?"available":"support");
            } else {
              print(images,0,"Select %s in image%s in interactive mode, from point (%g%s,%g%s,%g%s).",
                    feature_type==0?"point":feature_type==1?"segment":
                    feature_type==2?"rectangle":"ellipse",gmic_selection.data(),
                    value,sep=='%'?"%":"",value0,sep0=='%'?"%":"",value1,sep1=='%'?"%":"");

              unsigned int XYZ[3];
              cimg_forY(selection,l) {
                CImg<T> &img = images[selection[l]];
                XYZ[0] = (unsigned int)cimg::cut(cimg::round(sep=='%'?(img.width() - 1)*value/100:value),
                                                 0.,img.width() - 1.);
                XYZ[1] = (unsigned int)cimg::cut(cimg::round(sep0=='%'?(img.height() - 1)*value0/100:value0),
                                                 0.,img.height() - 1.);
                XYZ[2] = (unsigned int)cimg::cut(cimg::round(sep1=='%'?(img.depth() - 1)*value1/100:value1),
                                                 0.,img.depth() - 1.);
                if (gmic_display_window(0)) {
                  gmic_apply(select(gmic_display_window(0),feature_type,XYZ,
                                    (bool)exit_on_anykey,is_deep_selection),false);
                } else {
                  gmic_apply(select(images_names[selection[l]].data(),feature_type,XYZ,
                                    (bool)exit_on_anykey,is_deep_selection),false);
                }
              }
            }
          } else arg_error("select");
          is_change = true;
          ++position;
          continue;
        }

        // Serialize.
        if (!std::strcmp("serialize",command)) {
#define gmic_serialize(svalue_type,value_type) \
          if (!std::strcmp(argx,svalue_type)) \
            CImgList<value_type>(g_list,cimg::type<T>::string()==cimg::type<value_type>::string()). \
              get_serialize((bool)is_compressed).move_to(serialized);

          gmic_substitute_args(false);
#ifdef cimg_use_zlib
          bool is_compressed0 = true;
#else
          bool is_compressed0 = false;
#endif
          unsigned int is_compressed = is_compressed0?1U:0U, is_gmz = 1;
          if ((cimg_sscanf(argument,"%255[a-z123468 ]%c",
                           gmic_use_argx,&end)==1 ||
               cimg_sscanf(argument,"%255[a-z123468 ],%u%c",
                           argx,&is_compressed,&end)==2 ||
               cimg_sscanf(argument,"%255[a-z123468 ],%u,%u%c",
                           argx,&is_compressed,&is_gmz,&end)==3) &&
              (!std::strcmp(argx,"auto") ||
               !std::strcmp(argx,"uint8") || !std::strcmp(argx,"int8") ||
               !std::strcmp(argx,"uint16") || !std::strcmp(argx,"int16") ||
               !std::strcmp(argx,"uint32") || !std::strcmp(argx,"int32") ||
               !std::strcmp(argx,"uint64") || !std::strcmp(argx,"int64") ||
               !std::strcmp(argx,"float32") || !std::strcmp(argx,"float64")) &&
              is_compressed<=1 && is_gmz<=1) ++position;
          else { std::strcpy(argx,"auto"); is_compressed = is_compressed0?1U:0U; is_gmz = 1; }

          print(images,0,
                "Serialize %simage%s as %s%scompressed buffer%s, with datatype '%s'.",
                is_gmz?"named ":"",
                gmic_selection.data(),
                selection.height()>1?"":"a ",
                is_compressed?"":"un",
                selection.height()>1?"s":"",
                argx);

          if (selection) {
            g_list.assign(selection.height());
            CImgList<unsigned char> gmz_info;
            if (is_gmz) {
              gmz_info.assign(1 + selection.height());
              CImg<char>::string("GMZ").move_to(gmz_info[0]);
            }
            cimg_forY(selection,l) {
              uind = selection[l];
              g_list[l].assign(images[uind],true);
              if (is_gmz) CImg<char>::string(images_names[uind]).move_to(gmz_info[1 + l]);
            }
            if (is_gmz) (gmz_info>'x').unroll('y').move_to(g_list);
            if (!std::strcmp(argx,"auto")) std::strcpy(argx,CImg<T>::storage_type(g_list,false));
            CImg<unsigned char> serialized;
            gmic_serialize("uint8",cimg_uint8)
            else gmic_serialize("int8",cimg_int8)
              else gmic_serialize("uint16",cimg_uint16)
                else gmic_serialize("int16",cimg_int16)
                  else gmic_serialize("uint32",cimg_uint32)
                    else gmic_serialize("int32",cimg_int32)
                      else gmic_serialize("uint64",cimg_uint64)
                        else gmic_serialize("int64",cimg_int64)
                          else gmic_serialize("float32",cimg_float32)
                            else gmic_serialize("float64",cimg_float64)
                              else error(true,images,0,0,
                                         "Command 'serialize': Invalid "
                                         "specified pixel type '%s'.",
                                         argx);
            if (is_get) {
              serialized.move_to(images);
              images_names[selection[0]].get_copymark().move_to(images_names);
            } else {
              remove_images(images,images_names,selection,1,selection.height() - 1);
              serialized.move_to(images[selection[0]].assign());
            }
            g_list.assign();
          }
          is_change = true;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 't...'
        //-----------------------------
      gmic_commands_t :

        // Tangent.
        gmic_simple_command("tan",tan,"Compute pointwise tangent of image%s.");

        // Draw text.
        if (!std::strcmp("text",command)) {
          gmic_substitute_args(false);

          const char *p_argument = argument;
          is_cond = *argument==','; // is_cond = is_empty_string?
          if (is_cond) {
            CImg<char>((unsigned int)(std::strlen(argument) + 2),1,1,1).move_to(g_list_c);
            g_list_c(0,0) = ' ';
            std::strcpy(&g_list_c(0,1),argument);
            p_argument = g_list_c[0];
          }

          name.assign(4096);
          *argx = *argy = *argz = *name = *color = 0;
          double x = 0, y = 0, height = 16;
          bool is_custom_font = false;
          unsigned int nb_vals = 0;
          sep = sepx = sepy = sep0 = 0;
          opacity = 1;
          if ((cimg_sscanf(p_argument,"%4095[^,]%c",
                           name.data(),&end)==1 ||
               cimg_sscanf(p_argument,"%4095[^,],%255[0-9.eE%~+-]%c",
                           name.data(),gmic_use_argx,&end)==2 ||
               cimg_sscanf(p_argument,"%4095[^,],%255[0-9.eE%~+-],%255[0-9.eE%~+-]%c",
                           name.data(),argx,gmic_use_argy,&end)==3 ||
               cimg_sscanf(p_argument,"%4095[^,],%255[0-9.eE%~+-],%255[0-9.eE%~+-],%255[a-zA-Z_0-9.eE%+-]%c",
                           name.data(),argx,argy,gmic_use_argz,&end)==4 ||
               cimg_sscanf(p_argument,"%4095[^,],%255[0-9.eE%~+-],%255[0-9.eE%~+-],%255[a-zA-Z_0-9.eE%+-],%f%c",
                           name.data(),argx,argy,argz,&opacity,&end)==5 ||
               cimg_sscanf(p_argument,"%4095[^,],%255[0-9.eE%~+-],%255[0-9.eE%~+-],%255[a-zA-Z_0-9.eE%+-],%f,"
                           "%4095[0-9.eEinfa,+-]%c",
                           name.data(),argx,argy,argz,&opacity,gmic_use_color,&end)==6) &&
              (!*argx ||
               cimg_sscanf(argx,"%lf%c",&x,&end)==1 ||
               (cimg_sscanf(argx,"%lf%c%c",&x,&sepx,&end)==2 && (sepx=='%' || sepx=='~'))) &&
              (!*argy ||
               cimg_sscanf(argy,"%lf%c",&y,&end)==1 ||
               (cimg_sscanf(argy,"%lf%c%c",&y,&sepy,&end)==2 && (sepy=='%' || sepy=='~'))) &&
              (!*argz ||
               cimg_sscanf(argz,"%lf%c",&height,&end)==1 ||
               (cimg_sscanf(argz,"%lf%c%c",&height,&sep,&end)==2 && sep=='%') ||
               (is_custom_font = cimg::is_varname(argz))) &&
              height>=0) {

            if (!is_cond) {
              strreplace_fw(name);
              cimg::strunescape(name);
              if (*color) nb_vals = 1;
              for (const char *s = color; *s; ++s) if (*s==',') ++nb_vals;
            }

            if (is_custom_font) {
              print(images,0,"Draw text '%s' at position (%g%s,%g%s) on image%s, with font '%s', "
                    "opacity %g and color (%s).",
                    is_cond?"":name.data(),
                    x,sepx=='%'?"%":sepx=='~'?"~":"",
                    y,sepy=='%'?"%":sepy=='~'?"~":"",
                    gmic_selection.data(),
                    argz,opacity,
                    *color?color:"default");
              if (!is_cond) {
                CImgList<T> font;
                try {
                  unsigned int l_font = 0;
                  const CImg<char> s_font = get_variable(argz,variables_sizes,&images_names,&l_font);
                  CImgList<T>::get_unserialize(s_font,l_font + 1).move_to(font);
                } catch (CImgException&) {
                  error(true,images,0,0,
                        "Command 'text': Specified custom font '%s' is invalid.",
                        argz);
                }
                cimg_forY(selection,l) {
                  CImg<T> &img = images[selection[l]];
                  g_img.assign(std::max(img.spectrum(),(int)nb_vals),1,1,1,(T)0).fill_from_values(color,true);
                  gmic_apply(gmic_draw_text((float)x,(float)y,sepx,sepy,name,g_img,0,opacity,font,nb_vals),true);
                }
              }
            } else {
              print(images,0,"Draw text '%s' at position (%g%s,%g%s) on image%s, with font height %s, "
                    "opacity %g and color (%s).",
                    is_cond?"":name.data(),
                    x,sepx=='%'?"%":sepx=='~'?"~":"",
                    y,sepy=='%'?"%":sepy=='~'?"~":"",
                    gmic_selection.data(),
                    argz,opacity,
                    *color?color:"default");
              if (!is_cond)
                cimg_forY(selection,l) {
                  CImg<T> &img = images[selection[l]];
                  const unsigned int font_height = (unsigned int)cimg::round(sep=='%'?
                                                                             height*img.height()/100:height);
                  g_img.assign(std::max(img.spectrum(),(int)nb_vals),1,1,1,(T)0).fill_from_values(color,true);
                  gmic_apply(gmic_draw_text((float)x,(float)y,sepx,sepy,name,g_img,0,opacity,font_height,nb_vals),true);
                }
            }
          } else arg_error("text");
          g_img.assign();
          g_list_c.assign();
          is_change = true;
          ++position;
          continue;
        }

        // Tridiagonal solve.
        if (!std::strcmp("trisolve",command)) {
          gmic_substitute_args(true);
          sep = *indices = 0;
          if (cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep,&end)==2 &&
              sep==']' &&
              (ind=selection2cimg(indices,images.size(),images_names,"trisolve")).height()==1) {
            print(images,0,"Solve tridiagonal system AX = B, with B-vector%s and tridiagonal "
                  "A-matrix [%d].",
                  gmic_selection.data(),*ind);
            const CImg<double> A = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply_double(solve_tridiagonal(A));
          } else arg_error("trisolve");
          is_change = true;
          ++position;
          continue;
        }

        // Hyperbolic tangent.
        gmic_simple_command("tanh",tanh,"Compute pointwise hyperbolic tangent of image%s.");

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'u...'
        //-----------------------------
      gmic_commands_u :

        // Unroll.
        if (!std::strcmp("unroll",command)) {
          gmic_substitute_args(false);
          axis = 'y';
          if (is_xyzc(*argument) && !argument[1]) { axis = *argument; ++position; }
          print(images,0,"Unroll image%s along the '%c'-axis.",
                gmic_selection.data(),
                axis);
          cimg_forY(selection,l) gmic_apply(unroll(axis),false);
          is_change = true;
          continue;
        }

        // Remove custom command.
        if (!is_get && !std::strcmp("uncommand",item)) {
          gmic_substitute_args(false);
          if (*argument=='*' && !argument[1]) { // Discard all custom commands
            cimg::mutex(23);
            unsigned int nb_commands = 0;
            for (unsigned int i = 0; i<gmic_comslots; ++i) {
              nb_commands+=commands[i].size();
              commands[i].assign();
              commands_names[i].assign();
              commands_has_arguments[i].assign();
            }
            print(images,0,"Discard definitions of all custom commands (%u command%s).",
                  nb_commands,nb_commands>1?"s":"");
            cimg::mutex(23,0);
          } else { // Discard one or several custom command
            cimg::mutex(23);
            g_list_c = CImg<char>::string(argument).get_split(CImg<char>::vector(','),0,false);
            print(images,0,"Discard definition%s of custom command%s '%s'",
                  g_list_c.width()>1?"s":"",
                  g_list_c.width()>1?"s":"",
                  gmic_argument_text_printed());
            unsigned int nb_removed = 0;
            cimglist_for(g_list_c,l) {
              CImg<char>& arg_command = g_list_c[l];
              arg_command.resize(1,arg_command.height() + 1,1,1,0);
              strreplace_fw(arg_command);
              if (*arg_command) {
                hash = hashcode(arg_command,false);
                if (search_sorted(arg_command,commands_names[hash],commands_names[hash].size(),pattern)) {
                  commands_names[hash].remove(pattern);
                  commands[hash].remove(pattern);
                  commands_has_arguments[hash].remove(pattern);
                  ++nb_removed;
                }
              }
            }
            if (is_verbose) {
              cimg::mutex(29);
              unsigned int isiz = 0;
              for (unsigned int l = 0; l<gmic_comslots; ++l) isiz+=commands[l].size();
              std::fprintf(cimg::output()," (%u found, %u command%s left).",
                           nb_removed,isiz,isiz>1?"s":"");
              std::fflush(cimg::output());
              cimg::mutex(29,0);
            }
            g_list_c.assign();
            cimg::mutex(23,0);
          }
          ++position;
          continue;
        }

        // Unserialize.
        if (!std::strcmp("unserialize",command)) {
          print(images,0,"Unserialize image%s.",
                gmic_selection.data());
          int off = 0;
          cimg_forY(selection,l) {
            uind = selection[l] + off;
            const CImg<T>& img = gmic_check(images[uind]);
            CImgList<T>::get_unserialize(img).move_to(g_list);
            if (g_list) {
              const CImg<T>& back = g_list.back();
              if (back.width()==1 && back.depth()==1 && back.spectrum()==1 &&
                  back[0]=='G' && back[1]=='M' && back[2]=='Z' && !back[3]) { // .gmz serialization
                g_list_c.assign();
                const unsigned int pend = (unsigned int)back.size();
                for (unsigned int p = 4; p<pend; ) { // Retrieve list of image names
                  unsigned int np = p;
                  while (np<pend && back[np]) ++np;
                  if (np<pend) CImg<T>(back.data(p),1,++np - p,1,1,true).move_to(g_list_c);
                  p = np;
                }
                cimglist_for(g_list_c,q) g_list_c[q].unroll('x');
                if (g_list_c) g_list.remove();
              } else { // .cimg[z] serialization
                g_list_c.assign(g_list.size());
                if (is_get) images_names[uind].get_copymark().move_to(g_list_c[0]);
                else g_list_c[0] = images_names[uind];
                for (unsigned int i = 1; i<g_list_c.size(); ++i)
                  g_list_c[i - 1].get_copymark().move_to(g_list_c[i]);
              }
              if (g_list_c.width()>g_list.width())
                g_list_c.remove(g_list.width(),g_list_c.width() - 1);
              else if (g_list_c.width()<g_list.width())
                g_list_c.insert(g_list.width() - g_list_c.width(),CImg<char>::string("[unnamed]"));
              if (is_get) {
                g_list.move_to(images,~0U);
                g_list_c.move_to(images_names,~0U);
              } else {
                images.insert(g_list.size() - 1,uind + 1);
                images_names.insert(g_list.size() - 1,uind + 1);
                cimglist_for(g_list,i) {
                  g_list[i].move_to(images[uind + i]);
                  g_list_c[i].move_to(images_names[uind + i]);
                }
                off+=(int)g_list.size() - 1;
              }
            } else if (!is_get) {
              images.remove(uind);
              images_names.remove(uind);
              --off;
            }
          }
          g_list.assign();
          g_list_c.assign();
          is_change = true;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'v...'
        //-----------------------------
      gmic_commands_v :

        // Set verbosity
        // (actually only display a log message, since it has been already processed before).
        if (is_command_verbose) {
          if (*argument=='-' && !argument[1])
            print(images,0,"Decrement verbosity level (set to %d).",
                  verbosity);
          else if (*argument=='+' && !argument[1]) {
            if (is_very_verbose) print(images,0,"Increment verbosity level (set to %d).",
                                       verbosity);
          } else if ((verbosity>=1 && prev_verbosity>=1) || is_debug)
            print(images,0,"Set verbosity level to %d.",
                  verbosity);
          if (is_verbose_argument) ++position;
          continue;
        }

        // Vanvliet filter.
        if (!std::strcmp("vanvliet",command)) {
          gmic_substitute_args(false);
          unsigned int order = 0;
          float sigma = 0;
          axis = sep = 0;
          boundary = 1;
          if ((cimg_sscanf(argument,"%f,%u,%c%c",&sigma,&order,&axis,&end)==3 ||
               (cimg_sscanf(argument,"%f%c,%u,%c%c",&sigma,&sep,&order,&axis,&end)==4 &&
                sep=='%') ||
               cimg_sscanf(argument,"%f,%u,%c,%u%c",&sigma,&order,&axis,&boundary,&end)==4 ||
               (cimg_sscanf(argument,"%f%c,%u,%c,%u%c",
                            &sigma,&sep,&order,&axis,&boundary,&end)==5 && sep=='%')) &&
              sigma>=0 && order<=3 && is_xyzc(axis) && boundary<=3) {
            print(images,0,"Apply %u-order Vanvliet filter on image%s, along axis '%c' with standard "
                  "deviation %g%s and %s boundary conditions.",
                  order,gmic_selection.data(),axis,
                  sigma,sep=='%'?"%":"",
                  boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror");
            if (sep=='%') sigma = -sigma;
            cimg_forY(selection,l) gmic_apply(vanvliet(sigma,order,axis,boundary),true);
          } else arg_error("vanvliet");
          is_change = true;
          ++position;
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'w...'
        //-----------------------------
      gmic_commands_w :

        // While.
        if (!is_get && !std::strcmp("while",item)) {
          gmic_substitute_args(false);
          const CImg<char>& s = callstack.back();
          if (s[0]!='*' || s[1]!='d')
            error(true,images,0,0,
                  "Command 'while': Not associated to a 'do' command within the same scope.");
          is_cond = check_cond(argument,images,"while");
          if (is_very_verbose)
            print(images,0,"Reach 'while' command -> condition '%s' %s.",
                  gmic_argument_text_printed(),
                  is_cond?"holds":"does not hold");
          unsigned int *const dw = dowhiles.data(0,nb_dowhiles - 1);
          ++dw[1];
          if (is_cond) {
            position = dw[0] + 1;
            next_debug_line = dw[2];
            next_debug_filename = debug_filename;
            continue;
          } else {
            if (is_very_verbose) print(images,0,"End 'do...while' block.");
            --nb_dowhiles;
            callstack.remove();
            ++position;
            continue;
          }
        }

        // Warning.
        if (is_command_warn) {
          if (verbosity>=0 || is_debug || is_get) {
            gmic_substitute_args(false);
            bool force_visible = false;
            if ((*argument=='0' || *argument=='1') && argument[1]==',') {
              force_visible = *argument=='1'; argument+=2;
            }
            name.assign(argument,(unsigned int)std::strlen(argument) + 1);
            cimg::strunescape(name);
            const int _verbosity = ++verbosity;
            std::FILE *_file = 0;
            if (is_get) { _file = cimg::output(); verbosity = 1; cimg::output(stdout); }
            if (is_selection) warn(images,&selection,force_visible,"%s",name.data());
            else warn(images,&CImg<unsigned int>::empty(),force_visible,"%s",name.data());
            if (is_get) { verbosity = _verbosity; cimg::output(_file); }
            --verbosity;
          }
          ++position;
          continue;
        }

        // Display images in display window.
        sep = '0';
        if (!is_get &&
            (!std::strcmp("window",command) ||
             cimg_sscanf(command,"window%c%c",&sep,&end)==1 ||
             cimg_sscanf(command,"w%c%c",&sep,&end)==1) &&
            sep>='0' && sep<='9') {
          wind = (unsigned int)(sep - '0');
          gmic_substitute_args(false);
          int norm = -1, fullscreen = -1;
          float dimw = -1, dimh = -1, posx = cimg::type<float>::inf(), posy = cimg::type<float>::inf();
          sep0 = sep1 = sepx = sepy = *argx = *argy = *argz = *argc = *title = 0;
          if ((cimg_sscanf(argument,"%255[0-9.eE%+-]%c",
                           gmic_use_argx,&end)==1 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-]%c",
                           argx,gmic_use_argy,&end)==2 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%d%c",
                           argx,argy,&norm,&end)==3 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%d,%d%c",
                           argx,argy,&norm,&fullscreen,&end)==4 ||
               cimg_sscanf(argument,
                           "%255[0-9.eE%+-],%255[0-9.eE%+-],%d,%d,%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-]%c",
                           argx,argy,&norm,&fullscreen,gmic_use_argz,gmic_use_argc,&end)==6 ||
               cimg_sscanf(argument,
                           "%255[0-9.eE%+-],%255[0-9.eE%+-],%d,%d,%255[0-9.eE%+-],"
                           "%255[0-9.eE%+-],%255[^\n]",
                           argx,argy,&norm,&fullscreen,argz,argc,gmic_use_title)==7 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%d,%d,%255[^\n]",
                           &(*argx=*argz=*argc=0),argy,&norm,&fullscreen,title)==5 ||
               cimg_sscanf(argument,"%255[0-9.eE%+-],%255[0-9.eE%+-],%d,%255[^\n]",
                           argx,argy,&(norm=fullscreen=-1),title)==4 ||
               (norm=fullscreen=-1,cimg_sscanf(argument,
                                               "%255[0-9.eE%+-],%255[0-9.eE%+-],%255[^\n]",
                                               argx,argy,title))==3) &&
              (cimg_sscanf(argx,"%f%c",&dimw,&end)==1 ||
               (cimg_sscanf(argx,"%f%c%c",&dimw,&sep0,&end)==2 && sep0=='%')) &&
              (!*argy ||
               cimg_sscanf(argy,"%f%c",&dimh,&end)==1 ||
               (cimg_sscanf(argy,"%f%c%c",&dimh,&sep1,&end)==2 && sep1=='%')) &&
              (!*argz ||
               cimg_sscanf(argz,"%f%c",&posx,&end)==1 ||
               (cimg_sscanf(argz,"%f%c%c",&posx,&sepx,&end)==2 && sepx=='%')) &&
              (!*argc ||
               cimg_sscanf(argc,"%f%c",&posy,&end)==1 ||
               (cimg_sscanf(argc,"%f%c%c",&posy,&sepy,&end)==2 && sepy=='%')) &&
              (dimw>=0 || dimw==-1) &&
              (dimh>=0 || dimh==-1) &&
              norm>=-1 && norm<=3) ++position;
          else {
            dimw = dimh = -1;
            norm = fullscreen = -1;
            posx = posy = cimg::type<float>::inf();
            sep0 = sep1 = 0;
          }
          if (dimw==0 || dimh==0) dimw = dimh = 0;
          if (*title) { strreplace_fw(title); cimg::strunescape(title); }

          if (!is_display_available) {
            print(images,0,
                  "Display image%s in display window [%d] (skipped, no display %s).",
                  gmic_selection.data(),
                  wind,cimg_display?"available":"support");
          } else {

            // Get images to display and compute associated optimal size.
            unsigned int optw = 0, opth = 0;
            if (dimw && dimh) cimg_forY(selection,l) {
                const CImg<T>& img = gmic_check(images[selection[l]]);
                if (img) {
                  g_list.insert(img,~0U,true);
                  optw+=img._width + (img.depth()>1?img._depth:0U);
                  if (img.height()>(int)opth) opth = img._height + (img._depth>1?img._depth:0U);
                }
              }
            optw = optw?optw:sep0=='%'?CImgDisplay::screen_width():256;
            opth = opth?opth:sep1=='%'?CImgDisplay::screen_height():256;
            dimw = dimw<0?-1:cimg::round(sep0=='%'?optw*dimw/100:dimw);
            dimh = dimh<0?-1:cimg::round(sep1=='%'?opth*dimh/100:dimh);

            const bool is_move = !cimg::type<float>::is_inf(posx) && !cimg::type<float>::is_inf(posy);
            CImgDisplay &disp = gmic_display_window(wind);

            if (!dimw || !dimh) { // Close
              print(images,0,"Close display window [%d].",
                    wind);
              disp.assign();
            } else {
              if (disp) { // Update
                if (!selection) disp.show();
                disp.resize(dimw>0?(int)dimw:disp.window_width(),
                            dimh>0?(int)dimh:disp.window_height(),
                            false);
                if (is_move) {
                  if (sepx=='%') posx*=(CImgDisplay::screen_width() - disp.window_width())/100.f;
                  if (sepy=='%') posy*=(CImgDisplay::screen_height() - disp.window_height())/100.f;
                  disp.move((int)posx,(int)posy);
                }
                if (norm>=0) disp._normalization = (unsigned int)norm;
                if (*title && std::strcmp(disp.title(),title)) disp.set_title("%s",title);
                if (fullscreen>=0 && (bool)fullscreen!=disp.is_fullscreen()) disp.toggle_fullscreen(false);
              } else { // Create
                if (!*title) cimg_snprintf(title,_title.width(),"[G'MIC] Window #%u",wind);
                disp.assign(dimw>0?(int)dimw:optw,
                            dimh>0?(int)dimh:opth,
                            title,norm<0?3:norm,
                            fullscreen<0?false:(bool)fullscreen,
                            is_move);
                if (is_move) {
                  if (sepx=='%') posx*=(CImgDisplay::screen_width() - disp.window_width())/100.f;
                  if (sepy=='%') posy*=(CImgDisplay::screen_height() - disp.window_height())/100.f;
                  disp.move((int)posx,(int)posy);
                }
                if (norm==2) {
                  if (g_list) disp._max = (float)g_list.max_min(disp._min);
                  else { disp._min = 0; disp._max = 255; }
                }
              }
              if (is_move)
                print(images,0,
                      "Display image%s in %dx%d %sdisplay window [%d], "
                      "with%snormalization, "
                      "%sfullscreen, at position (%s,%s) and title '%s'.",
                      gmic_selection.data(),
                      disp.width(),disp.height(),disp.is_fullscreen()?"fullscreen ":"",
                      wind,
                      disp.normalization()==0?"out ":disp.normalization()==1?" ":
                      disp.normalization()==2?" 1st-time ":" auto-",
                      disp.is_fullscreen()?"":"no ",
                      argz,argc,
                      disp.title());
              else
                print(images,0,
                      "Display image%s in %dx%d %sdisplay window [%d], with%snormalization, "
                      "%sfullscreen and title '%s'.",
                      gmic_selection.data(),
                      disp.width(),disp.height(),disp.is_fullscreen()?"fullscreen ":"",
                      wind,
                      disp.normalization()==0?"out ":disp.normalization()==1?" ":
                      disp.normalization()==2?" 1st-time ":" auto-",
                      disp.is_fullscreen()?"":"no ",
                      disp.title());
              if (g_list) g_list.display(disp);
            }
            g_list.assign();
          }
          is_change = false;
          continue;
        }

        // Warp.
        if (!std::strcmp("warp",command)) {
          gmic_substitute_args(true);
          unsigned int mode = 0;
          double nb_frames = 1;
          interpolation = 1;
          boundary = 0;
          sep = 0;
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",
                            gmic_use_indices,&sep,&end)==2 && sep==']')||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u%c",
                           indices,&mode,&end)==2 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u%c",
                           indices,&mode,&interpolation,&end)==3 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%u%c",
                           indices,&mode,&interpolation,&boundary,&end)==4 ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u,%u,%u,%lf%c",
                           indices,&mode,&interpolation,&boundary,&nb_frames,&end)==5) &&
              (ind=selection2cimg(indices,images.size(),images_names,"warp")).height()==1 &&
              mode<=3 && interpolation<=2 && boundary<=3 && nb_frames>=0.5) {
            const CImg<T> warping_field = gmic_image_arg(*ind);
            nb_frames = cimg::round(nb_frames);
            if (nb_frames==1) {
              print(images,0,"Warp image%s with %s-%s displacement field [%u], %s interpolation, "
                    "%s boundary conditions.",
                    gmic_selection.data(),
                    mode<=2?"backward":"forward",(mode%2)?"relative":"absolute",
                    *ind,
                    interpolation==2?"cubic":interpolation==1?"linear":"nearest-neighbor",
                    boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror");
              cimg_forY(selection,l) gmic_apply(warp(warping_field,mode,interpolation,boundary),false);
            } else {
              print(images,0,"Warp image%s with %s-%s displacement field [%u], %s interpolation, "
                    "%s boundary conditions and %d frames.",
                    gmic_selection.data(),
                    mode<=2?"backward":"forward",(mode%2)?"relative":"absolute",
                    *ind,
                    interpolation==2?"cubic":interpolation==1?"linear":"nearest-neighbor",
                    boundary==0?"dirichlet":boundary==1?"neumann":boundary==2?"periodic":"mirror",
                    (int)nb_frames);
              unsigned int off = 0;
              CImg<T> _warp;
              if (!(mode%2)) _warp.assign(warping_field,false);

              cimg_forY(selection,l) {
                uind = selection[l] + off;
                CImg<T>& img = gmic_check(images[uind]);
                g_list.assign((int)nb_frames);
                cimglist_for(g_list,t)
                  if (mode%2) g_list[t] = img.get_warp(warping_field*(t/(nb_frames - 1)),mode,interpolation,boundary);
                  else {
                    cimg_forXYZ(_warp,x,y,z) {
                      const float fact = (float)(t/(nb_frames - 1));
                      if (_warp.spectrum()>0) _warp(x,y,z,0) = x + (warping_field(x,y,z,0) - x)*fact;
                      if (_warp.spectrum()>1) _warp(x,y,z,1) = y + (warping_field(x,y,z,1) - y)*fact;
                      if (_warp.spectrum()>2) _warp(x,y,z,2) = z + (warping_field(x,y,z,2) - z)*fact;
                    }
                    g_list[t] = img.get_warp(_warp,mode,interpolation,boundary);
                  }
                if (is_get) {
                  pattern = images_names.size();
                  images_names.insert((int)nb_frames);
                  images_names[uind].get_copymark().move_to(images_names[pattern]);
                  for (unsigned int i = 1; i<g_list.size(); ++i)
                    images_names[pattern + i - 1].get_copymark().move_to(images_names[pattern + i]);
                  g_list.move_to(images,~0U);
                } else {
                  images.insert((int)nb_frames - 1,uind + 1);
                  images_names.insert((int)nb_frames - 1,uind + 1);
                  g_list[0].move_to(images[uind]);
                  for (unsigned int i = 1; i<g_list.size(); ++i) {
                    g_list[i].move_to(images[uind + i]);
                    images_names[uind + i - 1].get_copymark().move_to(images_names[uind + i]);
                  }
                  off+=(int)nb_frames - 1;
                }
              }
              g_list.assign();
            }
          } else arg_error("warp");
          is_change = true;
          ++position;
          continue;
        }

        // Watershed transform.
        if (!std::strcmp("watershed",command)) {
          gmic_substitute_args(true);
          is_high_connectivity = 1;
          sep = *indices = 0;
          if (((cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sep,&end)==2 &&
                sep==']') ||
               cimg_sscanf(argument,"[%255[a-zA-Z0-9_.%+-]],%u%c",
                           indices,&is_high_connectivity,&end)==2) &&
              (ind=selection2cimg(indices,images.size(),images_names,"watershed")).height()==1 &&
              is_high_connectivity<=1) {
            print(images,0,"Compute watershed transform of image%s with priority map [%u] and "
                  "%s connectivity.",
                  gmic_selection.data(),*ind,is_high_connectivity?"high":"low");
            const CImg<T> priority = gmic_image_arg(*ind);
            cimg_forY(selection,l) gmic_apply(watershed(priority,(bool)is_high_connectivity),false);
          } else arg_error("watershed");
          is_change = true;
          ++position;
          continue;
        }

        // Wait for a given delay of for user events on display window.
        if (!is_get && !std::strcmp("wait",command)) {
          gmic_substitute_args(false);

          const char *const s_nodisplay = " (skipped, no display available).";
          double delay = 0;
          if (cimg_sscanf(argument,"%lf%c",&delay,&end)==1) ++position;
          else delay = 0;

          if (!is_selection) {
            if (is_display_available) { // Put all active windows in selection
              uind = 0;
              for (unsigned int l = 0; l<gmic_winslots; ++l) if (gmic_display_window(l)) ++uind;
              CImg<unsigned int>(1,uind).move_to(selection);
              uind = 0;
              for (unsigned int l = 0; l<gmic_winslots; ++l) if (gmic_display_window(l)) selection[uind++] = l;
            } else selection.assign();
            if (is_verbose) selection2string(selection,images_names,1,gmic_selection);
          }

          if (!delay) {
            print(images,0,"Wait for user events on display window%s%s.",
                  gmic_selection.data(),is_display_available?"":s_nodisplay);
            if (is_display_available) switch (selection.height()) {
              case 1 : CImgDisplay::wait(gmic_display_window(selection[0])); break;
              case 2 : CImgDisplay::wait(gmic_display_window(selection[0]),gmic_display_window(selection[1])); break;
              case 3 : CImgDisplay::wait(gmic_display_window(selection[0]),gmic_display_window(selection[1]),
                                         gmic_display_window(selection[2]));
                break;
              case 4 : CImgDisplay::wait(gmic_display_window(selection[0]),gmic_display_window(selection[1]),
                                         gmic_display_window(selection[2]),gmic_display_window(selection[3]));
                break;
              case 5 : CImgDisplay::wait(gmic_display_window(selection[0]),gmic_display_window(selection[1]),
                                         gmic_display_window(selection[2]),gmic_display_window(selection[3]),
                                         gmic_display_window(selection[4]));
                break;
              case 6 : CImgDisplay::wait(gmic_display_window(selection[0]),gmic_display_window(selection[1]),
                                         gmic_display_window(selection[2]),gmic_display_window(selection[3]),
                                         gmic_display_window(selection[4]),gmic_display_window(selection[5]));
                break;
              case 7 : CImgDisplay::wait(gmic_display_window(selection[0]),gmic_display_window(selection[1]),
                                         gmic_display_window(selection[2]),gmic_display_window(selection[3]),
                                         gmic_display_window(selection[4]),gmic_display_window(selection[5]),
                                         gmic_display_window(selection[6]));
                break;
              case 8 : CImgDisplay::wait(gmic_display_window(selection[0]),gmic_display_window(selection[1]),
                                         gmic_display_window(selection[2]),gmic_display_window(selection[3]),
                                         gmic_display_window(selection[4]),gmic_display_window(selection[5]),
                                         gmic_display_window(selection[6]),gmic_display_window(selection[7]));
                break;
              case 9 : CImgDisplay::wait(gmic_display_window(selection[0]),gmic_display_window(selection[1]),
                                         gmic_display_window(selection[2]),gmic_display_window(selection[3]),
                                         gmic_display_window(selection[4]),gmic_display_window(selection[5]),
                                         gmic_display_window(selection[6]),gmic_display_window(selection[7]),
                                         gmic_display_window(selection[8]));
                break;
              case 10 : CImgDisplay::wait(gmic_display_window(selection[0]),gmic_display_window(selection[1]),
                                          gmic_display_window(selection[2]),gmic_display_window(selection[3]),
                                          gmic_display_window(selection[4]),gmic_display_window(selection[5]),
                                          gmic_display_window(selection[6]),gmic_display_window(selection[7]),
                                          gmic_display_window(selection[8]),gmic_display_window(selection[9]));
                break;
              }

          } else if (delay<0) {
            delay = cimg::round(-delay);
            if (selection && is_display_available) {
              print(images,0,
                    "Flush display events of display window%s and wait for %g milliseconds%s.",
                    gmic_selection.data(),delay,is_display_available?"":s_nodisplay);
              is_cond = false;
              cimg_forY(selection,l) {
                CImgDisplay &disp = gmic_display_window(selection[l]);
                if (disp) {
                  disp.flush();
                  if (!is_cond) { disp.wait((unsigned int)delay); is_cond = true; }
                }
              }
            } else {
              print(images,0,
                    "Wait for %g milliseconds.",
                    delay);
              cimg::wait((unsigned int)delay);
            }

          } else { // delay>0
            delay = cimg::round(delay);
            if (selection && is_display_available) {
              print(images,0,"Wait for %g milliseconds, according to display window%s%s.",
                    delay,gmic_selection.data(),is_display_available?"":s_nodisplay);
              is_cond = false;
              cimg_forY(selection,l) {
                CImgDisplay &disp = gmic_display_window(selection[l]);
                if (disp && !is_cond) { disp.wait((unsigned int)delay); is_cond = true; }
              }
            } else {
              print(images,0,
                    "Sleep for %g milliseconds.",
                    delay);
              cimg::sleep((unsigned int)delay);
            }
          }
          continue;
        }

        goto gmic_commands_others;

        //-----------------------------
        // Commands starting by 'x...'
        //-----------------------------
      gmic_commands_x :

        // Bitwise xor.
        gmic_arithmetic_command("xor",
                                operator^=,
                                "Compute bitwise XOR of image%s by %g%s",
                                gmic_selection.data(),value,ssep,Tlong,
                                operator^=,
                                "Compute bitwise XOR of image%s by image [%d]",
                                gmic_selection.data(),ind[0],
                                operator_xoreq,
                                "Compute bitwise XOR of image%s by expression %s",
                                gmic_selection.data(),gmic_argument_text_printed(),
                                "Compute sequential bitwise XOR of image%s");

        goto gmic_commands_others;

        //----------------------------
        // Other (special) commands.
        //----------------------------
      gmic_commands_others :

        if (is_builtin_command) {

          // Left brace (not ignored previously, so starts a new generic code block).
          if (!is_get && *item=='{' && !item[1]) {
            if (is_debug_info && debug_line!=~0U) {
              gmic_use_argx;
              cimg_snprintf(argx,_argx.width(),"*block#%u",debug_line);
              CImg<char>::string(argx).move_to(callstack);
            } else CImg<char>::string("*block").move_to(callstack);
            continue;
          }

          // If...[elif]...[else]...endif.
          if (!is_get && (!std::strcmp("if",item) || (check_elif && !std::strcmp("elif",item)))) {
            gmic_substitute_args(false);
            is_cond = check_cond(argument,images,*item=='i'?"if":"elif");
            check_elif = false;
            if (*item=='i') {
              if (is_debug_info && debug_line!=~0U) {
                gmic_use_argx;
                cimg_snprintf(argx,_argx.width(),"*if#%u",debug_line);
                CImg<char>::string(argx).move_to(callstack);
              } else CImg<char>::string("*if").move_to(callstack);
              if (is_very_verbose) print(images,0,"Start 'if...endif' block -> condition '%s' %s.",
                                         gmic_argument_text_printed(),
                                         is_cond?"holds":"does not hold");
            } else if (is_very_verbose) print(images,0,"Reach 'elif' block -> condition '%s' %s.",
                                              gmic_argument_text_printed(),
                                              is_cond?"holds":"does not hold");
            if (!is_cond) {
              for (int nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
                it = commands_line[position];
                if (*it==1)
                  is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
                else {
                  it+=*it=='-';
                  if (!std::strcmp("if",it)) ++nb_levels;
                  else if (!std::strcmp("fi",it)) { if (!--nb_levels) --position; }
                  else if (nb_levels==1) {
                    if (!std::strcmp("else",it)) --nb_levels;
                    else if (!std::strcmp("elif",it)) { --nb_levels; check_elif = true; --position; }
                  }
                }
              }
              continue;
            }
            ++position;
            continue;
          }

          // Break or continue.
          bool is_continue = false;
          if (!is_get && (!std::strcmp("break",item) ||
                          (!std::strcmp("continue",item) && (is_continue=true)==true))) {
            const char
              *const com = is_continue?"continue":"break",
              *const Com = is_continue?"Continue":"Break";
            unsigned int callstack_repeat = 0, callstack_do = 0, callstack_for = 0, callstack_foreach = 0,
              callstack_local = 0;
            for (unsigned int l = callstack.size() - 1; l; --l) {
              const char *const s = callstack[l].data();
              const bool is_star = *s=='*';
              if (is_star && s[1]=='r') { callstack_repeat = l; break; }
              else if (is_star && s[1]=='d') { callstack_do = l; break; }
              else if (is_star && s[1]=='f') {
                if (s[4]!='e') callstack_for = l; else callstack_foreach = l;
                break;
              }
              else if (is_star && s[1]=='l') { callstack_local = l; break; }
              else if (!is_star || (s[1]!='i' && s[1]!='b')) break;
            }
            const char *stb = 0, *ste = 0;
            unsigned int callstack_ind = 0;
            int nb_levels = 0;
            if (callstack_repeat) {
              print(images,0,"%s %scurrent 'repeat...done' block.",
                    Com,is_continue?"to next iteration of ":"");
              for (nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
                it = commands_line[position];
                if (*it==1)
                  is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
                else {
                  _is_get = *it=='+';
                  it+=(_is_get || *it=='-');
                  gmic_if_flr ++nb_levels; gmic_elif_flr --nb_levels;
                }
              }
              callstack_ind = callstack_repeat;
              stb = "repeat"; ste = "done";
            } else if (callstack_do) {
              print(images,0,"%s %scurrent 'do...while' block.",
                    Com,is_continue?"to next iteration of ":"");
              for (nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
                it = commands_line[position];
                it+=*it=='-';
                if (!std::strcmp("do",it)) ++nb_levels;
                else if (!std::strcmp("while",it)) --nb_levels;
              }
              callstack_ind = callstack_do;
              stb = "do"; ste = "while";
            } else if (callstack_for) {
              print(images,0,"%s %scurrent 'for...done' block.",
                    Com,is_continue?"to next iteration of ":"");
              for (nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
                it = commands_line[position];
                if (*it==1)
                  is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
                else {
                  _is_get = *it=='+';
                  it+=(_is_get || *it=='-');
                  gmic_if_flr ++nb_levels; gmic_elif_flr --nb_levels;
                }
              }
              callstack_ind = callstack_for;
              stb = "for"; ste = "done";
            } else if (callstack_foreach) {
              print(images,0,"%s %scurrent 'foreach...done' block.",
                    Com,is_continue?"to next iteration of ":"");
              for (nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
                it = commands_line[position];
                if (*it==1)
                  is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
                else {
                  _is_get = *it=='+';
                  it+=(_is_get || *it=='-');
                  gmic_if_flr ++nb_levels; gmic_elif_flr --nb_levels;
                }
              }
              callstack_ind = callstack_foreach;
              stb = "foreach"; ste = "done";
            } else if (callstack_local) {
              print(images,0,"%s %scurrent local environment.",
                    Com,is_continue?"to end of ":"");
              for (nb_levels = 1; nb_levels && position<commands_line.size(); ++position) {
                it = commands_line[position];
                if (*it==1)
                  is_debug_info|=get_debug_info(commands_line[position].data(),next_debug_line,next_debug_filename);
                else {
                  _is_get = *it=='+';
                  it+=(_is_get || *it=='-');
                  gmic_if_flr ++nb_levels; gmic_elif_flr --nb_levels;
                }
              }
              callstack_ind = callstack_local;
              stb = "local"; ste = "done";
            } else {
              print(images,0,"%s",Com);
              error(true,images,0,0,
                    "Command '%s': There are no loops or local environment to %s.",com,com);
              continue;
            }

            if (nb_levels)
              error(true,images,0,0,
                    "Command '%s': Missing associated '%s' command.",stb,ste);
            if (is_continue || callstack_local || callstack_foreach) {
              if (callstack_foreach && !is_continue) { // Break 'foreach...done' loop
                unsigned int *const fed = foreachdones.data(0,nb_foreachdones - 1);
                fed[0]+=fed[1]; // Force loop to end at next 'done'
                fed[1] = 0;
              }
              if (callstack_ind<callstack.size() - 1) callstack.remove(callstack_ind + 1,callstack.size() - 1);
              --position;
            } else {
              callstack.remove(callstack_ind,callstack.size() - 1);
              if (callstack_do) { --nb_dowhiles; ++position; }
              else if (callstack_repeat) --nb_repeatdones;
              else --nb_fordones;
            }
            continue;
          }

          // Compute direct or inverse FFT.
          const bool inv_fft = !std::strcmp("ifft",command);
          if (!std::strcmp("fft",command) || inv_fft) {
            gmic_substitute_args(false);
            bool is_valid_argument = *argument!=0;
            if (is_valid_argument) for (const char *s = argument; *s; ++s) {
                const char _s = *s;
                if (_s!='x' && _s!='y' && _s!='z') { is_valid_argument = false; break; }
              }
            if (is_valid_argument) {
              print(images,0,"Compute %sfourier transform of image%s along the '%s'-ax%cs with complex pair%s",
                    inv_fft?"inverse ":"",
                    gmic_selection.data(),
                    gmic_argument_text_printed(),
                    std::strlen(argument)>1?'e':'i',
                    selection.height()>2?"s":selection.height()>=1?"":"().");
              ++position;
            } else
              print(images,0,"Compute %sfourier transform of image%s with complex pair%s",
                    inv_fft?"inverse ":"",
                    gmic_selection.data(),
                    selection.height()>2?"s":selection.height()>=1?"":" ().");
            cimg_forY(selection,l) {
              const unsigned int
                uind0 = selection[l],
                uind1 = l + 1<selection.height()?selection[l + 1]:~0U;
              CImg<T> &img0 = gmic_check(images[uind0]),
                      &img1 = uind1!=~0U?gmic_check(images[uind1]):CImg<T>::empty();
              if (uind1!=~0U) { // Complex transform
                if (is_verbose) {
                  cimg::mutex(29);
                  std::fprintf(cimg::output()," ([%u],[%u])%c",uind0,uind1,
                               l>=selection.height() - 2?'.':',');
                  std::fflush(cimg::output());
                  cimg::mutex(29,0);
                }
                if (is_get) {
                  g_list.assign(img0,img1);
                  if (is_valid_argument) for (const char *s = argument; *s; ++s) g_list.FFT(*s,inv_fft);
                  else g_list.FFT(inv_fft);
                  g_list.move_to(images,~0U);
                  images_names[uind0].get_copymark().move_to(images_names);
                  images_names.back().get_copymark().move_to(images_names);
                } else {
                  g_list.assign(2);
                  g_list[0].swap(img0);
                  g_list[1].swap(img1);
                  if (is_valid_argument) for (const char *s = argument; *s; ++s) g_list.FFT(*s,inv_fft);
                  else g_list.FFT(inv_fft);
                  g_list[0].swap(img0);
                  g_list[1].swap(img1);
                  images_names[uind0].get_copymark().move_to(images_names[uind1]);
                }
                ++l;
              } else { // Real transform
                if (is_verbose) {
                  cimg::mutex(29);
                  std::fprintf(cimg::output()," ([%u],0)%c",uind0,
                               l>=selection.height() - 2?'.':',');
                  std::fflush(cimg::output());
                  cimg::mutex(29,0);
                }
                if (is_get) {
                  g_list.assign(img0);
                  CImg<T>(g_list[0].width(),g_list[0].height(),g_list[0].depth(),g_list[0].spectrum(),(T)0).
                    move_to(g_list);
                  if (is_valid_argument) for (const char *s = argument; *s; ++s) g_list.FFT(*s,inv_fft);
                  else g_list.FFT(inv_fft);
                  g_list.move_to(images,~0U);
                  images_names[uind0].get_copymark().move_to(images_names);
                  images_names.back().get_copymark().move_to(images_names);
                } else {
                  g_list.assign(1);
                  g_list[0].swap(img0);
                  CImg<T>(g_list[0].width(),g_list[0].height(),g_list[0].depth(),g_list[0].spectrum(),(T)0).
                    move_to(g_list);
                  if (is_valid_argument) for (const char *s = argument; *s; ++s) g_list.FFT(*s,inv_fft);
                  else g_list.FFT(inv_fft);
                  g_list[0].swap(img0);
                  g_list[1].move_to(images,uind0 + 1);
                  images_names[uind0].get_copymark().move_to(images_names,uind0 + 1);
                }
              }
            }
            g_list.assign();
            is_change = true;
            continue;
          }

          // Rescale a 3D object (* or /).
          const bool divide3d = !std::strcmp("div3d",command);
          if (!std::strcmp("mul3d",command) || divide3d) {
            gmic_substitute_args(false);
            float sx = 0, sy = 1, sz = 1;
            if ((cimg_sscanf(argument,"%f%c",
                             &sx,&end)==1 && ((sz=sy=sx),1)) ||
                cimg_sscanf(argument,"%f,%f%c",
                            &sx,&sy,&end)==2 ||
                cimg_sscanf(argument,"%f,%f,%f%c",
                            &sx,&sy,&sz,&end)==3) {
              if (divide3d)
                print(images,0,"Scale 3D object%s with factors (1/%g,1/%g,1/%g).",
                      gmic_selection.data(),
                      sx,sy,sz);
              else
                print(images,0,"Scale 3D object%s with factors (%g,%g,%g).",
                      gmic_selection.data(),
                      sx,sy,sz);
              cimg_forY(selection,l) {
                uind = selection[l];
                CImg<T>& img = images[uind];
                try {
                  if (divide3d) { gmic_apply(scale_CImg3d(1/sx,1/sy,1/sz),true); }
                  else gmic_apply(scale_CImg3d(sx,sy,sz),true);
                } catch (CImgException&) {
                  if (!img.is_CImg3d(true,&(*gmic_use_message=0)))
                    error(true,images,0,0,
                          "Command '%s3d': Invalid 3D object [%d], in image%s (%s).",
                          divide3d?"div":"mul",uind,gmic_selection_err.data(),message);
                  else throw;
                }
              }
            } else { if (divide3d) arg_error("div3d"); else arg_error("mul3d"); }
            is_change = true;
            ++position;
            continue;
          }
        } // if (is_builtin_command)

        // Execute custom command.
        if (!is_command_input && is_command) {
          if (hash_custom==~0U) { // Probably a '+builtin_command' not supporting '+' prepend (e.g. +v)
            hash_custom = hashcode(command,false);
            is_command = search_sorted(command,commands_names[hash_custom],
                                       commands_names[hash_custom].size(),ind_custom);
          }

          if (is_command) {
            bool has_arguments = false, _is_noarg = false;
            CImg<char> substituted_command(1024);
            char *ptr_sub = substituted_command.data();
            const char
              *const command_code = commands[hash_custom][ind_custom].data(),
              *const command_code_back = &commands[hash_custom][ind_custom].back(),
              *const command_name = is_specialized_get?_command:command;

            if (is_debug) {
              CImg<char> command_code_text(264);
              const unsigned int ls = (unsigned int)std::strlen(command_code);
              if (ls>=264) {
                std::memcpy(command_code_text.data(),command_code,128);
                std::memcpy(command_code_text.data() + 128,"(...)",5);
                std::memcpy(command_code_text.data() + 133,command_code + ls - 130,131);
              } else std::strcpy(command_code_text.data(),command_code);
              for (char *ptrs = command_code_text, *ptrd = ptrs; *ptrs || (bool)(*ptrd=0);
                   ++ptrs)
                if (*ptrs==1) do ++ptrs; while (*ptrs!=' '); else *(ptrd++) = *ptrs;
              debug(images,"Found custom command '%s: %s' (%s).",
                    command_name,command_code_text.data(),
                    commands_has_arguments[hash_custom](ind_custom,0)?"takes arguments":
                    "takes no arguments");
            }

            CImgList<char> arguments(32);
            // Set $0 to be the command name.
            CImg<char>::string(command).move_to(arguments[0]);
            unsigned int nb_arguments = 0;

            if (commands_has_arguments[hash_custom](ind_custom,0)) { // Command takes arguments
              gmic_substitute_args(false);

              // Extract possible command arguments.
              for (const char *ss = argument, *_ss = ss; _ss; ss =_ss + 1)
                if ((_ss=std::strchr(ss,','))!=0) {
                  if (ss==_ss) ++nb_arguments;
                  else {
                    if (++nb_arguments>=arguments.size())
                      arguments.insert(2 + 2*nb_arguments - arguments.size());
                    CImg<char> arg_item(ss,(unsigned int)(_ss - ss + 1));
                    arg_item.back() = 0;
                    arg_item.move_to(arguments[nb_arguments]);
                  }
                } else {
                  if (*ss) {
                    if (++nb_arguments>=arguments.size())
                      arguments.insert(1 + nb_arguments - arguments.size());
                    if (*ss!=',') CImg<char>::string(ss).move_to(arguments[nb_arguments]);
                  }
                  break;
                }

              if (is_debug) {
                debug(images,"Found %d given argument%s for command '%s'%s",
                      nb_arguments,nb_arguments!=1?"s":"",
                      command_name,nb_arguments>0?":":".");
                for (unsigned int i = 1; i<=nb_arguments; ++i)
                  if (arguments[i]) debug(images,"  $%d = '%s'",i,arguments[i].data());
                  else debug(images,"  $%d = (undefined)",i);
              }
            }

            // Substitute arguments in custom command expression.
            CImg<char> inbraces;

            for (const char *nsource = command_code; *nsource;)
              if (*nsource!='$') {

                // If not starting with '$'.
                const char *const nsource0 = nsource;
                nsource = std::strchr(nsource0,'$');
                if (!nsource) nsource = command_code_back;
                CImg<char>(nsource0,(unsigned int)(nsource - nsource0),1,1,1,true).
                  append_string_to(substituted_command,ptr_sub);
              } else { // '$' expression found
                CImg<char> substr(324);
                inbraces.assign(1,1,1,1,0);
                int iind = 0, iind1 = 0, l_inbraces = 0;
                bool is_braces = false;
                sep = 0;

                if (nsource[1]=='{') {
                  const char *const ptr_beg = nsource + 2, *ptr_end = ptr_beg;
                  unsigned int p = 0;
                  for (p = 1; p>0 && *ptr_end; ++ptr_end) {
                    if (*ptr_end=='{') ++p;
                    if (*ptr_end=='}') --p;
                  }
                  if (p) { CImg<char>::append_string_to(*(nsource++),substituted_command,ptr_sub); continue; }
                  l_inbraces = (int)(ptr_end - ptr_beg - 1);
                  if (l_inbraces>0) inbraces.assign(ptr_beg,l_inbraces + 1).back() = 0;
                  is_braces = true;
                }

                // Substitute $# -> maximum index of known arguments.
                if (nsource[1]=='#') {
                  nsource+=2;
                  cimg_snprintf(substr,substr.width(),"%u",nb_arguments);
                  CImg<char>(substr.data(),(unsigned int)std::strlen(substr),1,1,1,true).
                    append_string_to(substituted_command,ptr_sub);
                  has_arguments = true;

                  // Substitute $* -> copy of the specified arguments string.
                } else if (nsource[1]=='*') {
                  nsource+=2;
                  CImg<char>(argument,(unsigned int)std::strlen(argument),1,1,1,true).
                    append_string_to(substituted_command,ptr_sub);
                  has_arguments = true;

                  // Substitute $"*" -> copy of the specified "quoted" arguments string.
                } else if (nsource[1]=='\"' && nsource[2]=='*' && nsource[3]=='\"') {
                  nsource+=4;
                  for (unsigned int i = 1; i<=nb_arguments; ++i) {
                    CImg<char>::append_string_to('\"',substituted_command,ptr_sub);
                    CImg<char>(arguments[i].data(),arguments[i].width() - 1,1,1,1,true).
                      append_string_to(substituted_command,ptr_sub);
                    CImg<char>::append_string_to('\"',substituted_command,ptr_sub);
                    if (i!=nb_arguments) CImg<char>::append_string_to(',',substituted_command,ptr_sub);
                  }
                  has_arguments = true;

                  // Substitute $[] -> List of selected image indices.
                } else if (nsource[1]=='[' && nsource[2]==']') {
                  nsource+=3;
                  cimg_forY(selection,i) {
                    cimg_snprintf(substr,substr.width(),"%u,",selection[i]);
                    CImg<char>(substr.data(),(unsigned int)std::strlen(substr),1,1,1,true).
                      append_string_to(substituted_command,ptr_sub);
                  }
                  if (selection) --ptr_sub;

                  // Substitute $= -> transfer (quoted) arguments to named variables.
                } else if (nsource[1]=='=' &&
                           cimg_sscanf(nsource + 2,"%255[a-zA-Z0-9_]",gmic_use_title)==1 &&
                           (*title<'0' || *title>'9')) {
                  nsource+=2 + std::strlen(title);
                  CImg<char>::append_string_to(' ',substituted_command,ptr_sub);
                  for (unsigned int i = 0; i<=nb_arguments; ++i) {
                    cimg_snprintf(substr,substr.width(),"%s%u%c",title,i,i==nb_arguments?'=':',');
                    CImg<char>(substr.data(),(unsigned int)std::strlen(substr),1,1,1,true).
                      append_string_to(substituted_command,ptr_sub);
                  }
                  for (unsigned int i = 0; i<=nb_arguments; ++i) {
                    CImg<char>::append_string_to('\"',substituted_command,ptr_sub);
                    CImg<char>(arguments[i].data(),arguments[i].width() - 1,1,1,1,true).
                      append_string_to(substituted_command,ptr_sub);
                    CImg<char>::append_string_to('\"',substituted_command,ptr_sub);
                    CImg<char>::append_string_to(i==nb_arguments?' ':',',substituted_command,ptr_sub);
                  }
                  has_arguments = true;

                  // Substitute $i and ${i} -> value of the i^th argument.
                } else if ((cimg_sscanf(nsource,"$%d",&iind)==1 ||
                            (cimg_sscanf(nsource,"${%d%c",&iind,&sep)==2 && sep=='}'))) {
                  const int niind = iind + (iind<0?(int)nb_arguments + 1:0);
                  if ((niind<=0 && iind) || niind>=arguments.width() || !arguments[niind]) {
                    error(true,images,0,command_name,
                          "Command '%s': Undefined argument '$%d', in expression '$%s%d%s' "
                          "(for %u argument%s specified).",
                          command_name,iind,sep=='}'?"{":"",iind,sep=='}'?"}":"",
                          nb_arguments,nb_arguments!=1?"s":"");
                  }
                  nsource+=cimg_snprintf(substr,substr.width(),"$%d",iind) + (sep=='}'?2:0);
                  if (arguments[niind].width()>1)
                    CImg<char>(arguments[niind].data(),arguments[niind].width() - 1,1,1,1,true).
                      append_string_to(substituted_command,ptr_sub);
                  if (niind!=0) has_arguments = true;

                  // Substitute ${i=$j} -> value of the i^th argument, or the default value,
                  // i.e. the value of another argument.
                } else if (cimg_sscanf(nsource,"${%d=$%d%c",&iind,&iind1,&sep)==3 && sep=='}' &&
                           iind>0) {
                  const int niind1 = iind1 + (iind1<0?(int)nb_arguments + 1:0);
                  if (niind1<=0 || niind1>=arguments.width() || !arguments[niind1])
                    error(true,images,0,command_name,
                          "Command '%s': Undefined argument '$%d', in expression '${%d=$%d}' "
                          "(for %u argument%s specified).",
                          command_name,iind1,iind,iind1,
                          nb_arguments,nb_arguments!=1?"s":"");
                  nsource+=cimg_snprintf(substr,substr.width(),"${%d=$%d}",iind,iind1);
                  if (iind>=arguments.width()) arguments.insert(2 + 2*iind - arguments.size());
                  if (!arguments[iind]) {
                    arguments[iind] = arguments[niind1];
                    if (iind>(int)nb_arguments) nb_arguments = (unsigned int)iind;
                  }
                  if (arguments[iind].width()>1)
                    CImg<char>(arguments[iind].data(),arguments[iind].width() - 1,1,1,1,true).
                      append_string_to(substituted_command,ptr_sub);
                  has_arguments = true;

                  // Substitute ${i=$#} -> value of the i^th argument, or the default value,
                  // i.e. the maximum index of known arguments.
                } else if (cimg_sscanf(nsource,"${%d=$#%c",&iind,&sep)==2 && sep=='}' &&
                           iind>0) {
                  if (iind>=arguments.width()) arguments.insert(2 + 2*iind - arguments.size());
                  if (!arguments[iind]) {
                    cimg_snprintf(substr,substr.width(),"%u",nb_arguments);
                    CImg<char>::string(substr).move_to(arguments[iind]);
                    if (iind>(int)nb_arguments) nb_arguments = (unsigned int)iind;
                  }
                  nsource+=cimg_snprintf(substr,substr.width(),"${%d=$#}",iind);
                  if (arguments[iind].width()>1)
                    CImg<char>(arguments[iind].data(),arguments[iind].width() - 1,1,1,1,true).
                      append_string_to(substituted_command,ptr_sub);
                  has_arguments = true;

                  // Substitute ${i=default} -> value of the i^th argument,
                  // or the specified default value.
                } else if (cimg_sscanf(inbraces,"%d%c",&iind,&sep)==2 && sep=='=' &&
                           iind>0) {
                  nsource+=l_inbraces + 3;
                  if (iind>=arguments.width()) arguments.insert(2 + 2*iind - arguments.size());
                  if (!arguments[iind]) {
                    CImg<char>::string(inbraces.data() +
                                       cimg_snprintf(substr,substr.width(),"%d=",iind)).
                      move_to(arguments[iind]);
                    if (iind>(int)nb_arguments) nb_arguments = (unsigned int)iind;
                  }
                  if (arguments[iind].width()>1)
                    CImg<char>(arguments[iind].data(),arguments[iind].width() - 1,1,1,1,true).
                      append_string_to(substituted_command,ptr_sub);
                  has_arguments = true;

                  // Substitute any other expression starting by '$'.
                } else {

                  // Substitute ${subset} -> values of the selected subset of arguments,
                  // separated by ','.
                  bool is_valid_subset = false;
                  if (is_braces) {
                    const char c = *inbraces, nc = c?inbraces[1]:0;
                    if (c=='^' || c==':' || c=='.' || (c>='0' && c<='9') ||
                        (c=='-' && !((nc>='a' && nc<='z') ||
                                     (nc>='A' && nc<='Z') ||
                                     nc=='_'))) {

                      CImg<unsigned int> inds;
                      status.move_to(o_status); // Save status because 'selection2cimg' can change it
                      const int o_verbosity = verbosity;
                      const bool o_is_debug = is_debug;
                      verbosity = 0;
                      is_debug = false;
                      try {
                        inds = selection2cimg(inbraces,nb_arguments + 1,
                                              CImgList<char>::empty(),"",false);
                        is_valid_subset = true;
                      } catch (...) { inds.assign(); is_valid_subset = false; }
                      is_debug = o_is_debug;
                      verbosity = o_verbosity;
                      o_status.move_to(status);

                      if (is_valid_subset) {
                        nsource+=l_inbraces + 3;
                        if (inds) {
                          cimg_forY(inds,j) {
                            uind = inds[j];
                            if (uind) has_arguments = true;
                            if (!arguments[uind])
                              error(true,images,0,command_name,
                                    "Command '%s': Undefined argument '$%d', "
                                    "in expression '${%s}'.",
                                    command_name,uind,inbraces.data());
                            CImg<char>(arguments[uind],true).append_string_to(substituted_command,ptr_sub);
                            *(ptr_sub - 1) = ',';
                          }
                          --ptr_sub;
                          has_arguments = true;
                        }
                      }
                    }
                  }
                  if (!is_valid_subset) CImg<char>::append_string_to(*(nsource++),substituted_command,ptr_sub);
                }
              }
            *ptr_sub = 0;

            // Substitute special character codes appearing outside strings.
            bool is_dquoted = false, is_escaped = false;
            for (char *s = substituted_command.data(); *s; ++s) {
              const char c = *s;
              if (is_escaped) is_escaped = false;
              else if (c=='\\') is_escaped = true;
              else if (c=='\"') is_dquoted = !is_dquoted;
              else if (!is_dquoted) _strreplace_fw(*s);
            }

            if (is_debug) {
              CImg<char> command_code_text(264);
              const unsigned int l = (unsigned int)std::strlen(substituted_command.data());
              if (l>=264) {
                std::memcpy(command_code_text.data(),substituted_command.data(),128);
                std::memcpy(command_code_text.data() + 128,"(...)",5);
                std::memcpy(command_code_text.data() + 133,substituted_command.data() + l - 130,131);
              } else std::strcpy(command_code_text.data(),substituted_command.data());
              for (char *ptrs = command_code_text, *ptrd = ptrs; *ptrs || (bool)(*ptrd=0);
                   ++ptrs)
                if (*ptrs==1) do ++ptrs; while (*ptrs!=' '); else *(ptrd++) = *ptrs;
              debug(images,"Expand command line for command '%s' to: '%s'.",
                    command_name,command_code_text.data());
            }

            const CImgList<char>
              ncommands_line = commands_line_to_CImgList(substituted_command.data());
            CImg<unsigned int> nvariables_sizes(gmic_varslots);
            cimg_forX(nvariables_sizes,l) nvariables_sizes[l] = variables[l]->size();
            g_list.assign(selection.height());
            g_list_c.assign(selection.height());

            unsigned int nposition = 0;
            gmic_exception exception;
            const unsigned int
              previous_debug_filename = debug_filename,
              previous_debug_line = debug_line;
            CImg<char>::string(commands_names[hash_custom][ind_custom]).move_to(callstack);

            if (is_get) { // Call to '+command'
              cimg_forY(selection,l) {
                uind = selection[l];
                g_list[l] = images[uind];
                g_list_c[l] = images_names[uind];
              }

              try {
                is_debug_info = false;
                --verbosity;
                _run(ncommands_line,nposition,g_list,g_list_c,images,images_names,nvariables_sizes,&_is_noarg,
                     argument,&selection,false);
                ++verbosity;
              } catch (gmic_exception &e) {
                cimg::swap(exception._command,e._command);
                cimg::swap(exception._message,e._message);
              }

              g_list.move_to(images,~0U);
              cimglist_for(g_list_c,l) g_list_c[l].get_copymark().move_to(g_list_c[l]);

              g_list_c.move_to(images_names,~0U);
            } else { // Call to 'command'
              cimg::mutex(27);
              cimg_forY(selection,l) {
                uind = selection[l];
                if ((images[uind].width() || images[uind].height()) && !images[uind].spectrum()) {
                  selection2string(selection,images_names,1,name);
                  error(true,images,0,0,
                        "Command '%s': Invalid selection%s "
                        "(image [%u] is already used in another thread).",
                        command_name,name.data() + (*name=='s'?1:0),uind);
                }
                if (images[uind].is_shared())
                  g_list[l].assign(images[uind],false);
                else {
                  g_list[l].swap(images[uind]);
                  // Small hack to be able to track images of the selection passed to the new environment.
                  std::memcpy(&images[uind]._width,&g_list[l]._data,sizeof(void*));
                  images[uind]._spectrum = 0;
                }
                g_list_c[l] = images_names[uind]; // Make a copy to be still able to recognize 'pass[label]'
              }
              cimg::mutex(27,0);

              try {
                is_debug_info = false;
                --verbosity;
                _run(ncommands_line,nposition,g_list,g_list_c,images,images_names,nvariables_sizes,&_is_noarg,
                     argument,&selection,false);
                ++verbosity;
              } catch (gmic_exception &e) {
                cimg::swap(exception._command,e._command);
                cimg::swap(exception._message,e._message);
              }
              if (run_main_) { --verbosity; run_main_ = false; }

              const unsigned int nb = std::min((unsigned int)selection.height(),g_list.size());
              if (nb>0) {
                for (unsigned int i = 0; i<nb; ++i) {
                  uind = selection[i];
                  if (images[uind].is_shared()) {
                    images[uind] = g_list[i];
                    g_list[i].assign();
                  } else images[uind].swap(g_list[i]);
                  images_names[uind].swap(g_list_c[i]);
                }
                g_list.remove(0,nb - 1);
                g_list_c.remove(0,nb - 1);
              }
              if (nb<(unsigned int)selection.height())
                remove_images(images,images_names,selection,nb,selection.height() - 1);
              else if (g_list) {
                const unsigned uind0 = selection && !is_specialized_get?selection.back() + 1:images.size();
                g_list_c.move_to(images_names,uind0);
                g_list.move_to(images,uind0);
              }
            }
            for (unsigned int l = 0; l<gmic_varslots/2; ++l) if (variables[l]->size()>nvariables_sizes[l]) {
                if (variables_lengths[l]->_width - nvariables_sizes[l]>variables_lengths[l]->_width/2)
                  variables_lengths[l]->resize(nvariables_sizes[l],1,1,1,0);
                variables_names[l]->remove(nvariables_sizes[l],variables[l]->size() - 1);
                variables[l]->remove(nvariables_sizes[l],variables[l]->size() - 1);
              }
            callstack.remove();
            debug_filename = previous_debug_filename;
            debug_line = previous_debug_line;
            is_return = false;
            g_list.assign();
            g_list_c.assign();
            if (has_arguments && !_is_noarg) ++position;
            if (exception._message) throw exception;
            continue;
          }
        }
      } // if (is_command) {

      // Variable assignment.
      if (!is_command_input && (*item=='_' || (*item>='a' && *item<='z') || (*item>='A' && *item<='Z'))) {
        const char *const s_equal = std::strchr(item + 1,'=');
        if (s_equal) {
          const char *s_end_left = s_equal;
          sep0 = *(s_equal - 1);
          sep1 = s_equal>item + 1?*(s_equal - 2):0;
          if (sep1=='.' && sep0==sep1) {
            s_end_left = s_equal - 2;
            sep0 = ',';
          } else if ((sep1=='>' || sep1=='<') && sep0==sep1)
            s_end_left = s_equal - 2;
          else if (sep0=='+' || sep0=='-' || sep0=='*' || sep0=='/' || sep0=='.' ||
                   sep0=='%' || sep0=='&' || sep0=='|' || sep0=='^' || sep0==':')
            s_end_left = s_equal - 1;
          else
            sep0 = '=';
          is_cond = sep0=='+' || sep0=='-' || sep0=='*' || sep0=='/' || sep0=='%' || sep0=='&' || sep0=='|' ||
            sep0=='^' || sep0==':' || sep0=='<' || sep0=='>'; // Right-hand side must be evaluated as a math expression?
          const char *const s_operation = sep0=='='?"":sep0==':'?":":sep0=='+'?"+":sep0=='-'?"-":sep0=='*'?"*":
            sep0=='/'?"/":sep0=='%'?"%":sep0=='&'?"&":sep0=='|'?"|":sep0=='^'?"^":sep0=='<'?"<<":">>";

          sep = s_end_left>item?*(s_end_left - 1):0;
          if ((sep>='a' && sep<='z') || (sep>='A' && sep<='Z') || (sep>='0' && sep<='9') || sep=='_') {

            // Check validity of variable name(s).
            CImgList<char> varnames, varvalues;
            bool is_valid_expr = true;
            const char *s = std::strchr(item,',');

            if (!s || s>s_equal) { // Single variable name
              is_valid_expr = cimg_sscanf(item,"%255[a-zA-Z0-9_]",gmic_use_title)==1 && (*title<'0' || *title>'9');
              is_valid_expr&=(item + std::strlen(title)==s_end_left);
              if (is_valid_expr) {
                s = s_equal + 1;
                varnames.insert(1).back().assign(title,(unsigned int)std::strlen(title) + 1,1,1,1,true);
                varvalues.insert(1).back().assign(s,(unsigned int)std::strlen(s) + 1,1,1,1,true);
              }
            } else { // Multiple variable names

              // Split list of variable names in 'varnames'.
              gmic_use_title;
              s = item;
              while (s<s_end_left) {
                *title = 0;
                const char *ns = std::strchr(s,',');
                if (ns==s) { is_valid_expr = false; break; }
                if (!ns || ns>s_equal) ns = s_end_left;
                CImg<char>(s,(unsigned int)(ns - s + 1)).move_to(name);
                name.back() = 0;
                if (cimg_sscanf(name,"%255[a-zA-Z0-9_]%c",title,&sep)==1 && (*title<'0' || *title>'9'))
                  name.move_to(varnames);
                else { is_valid_expr = false; break; }
                s = ns + 1;
              }
              if (!is_valid_expr) { // Variable names are invalid
                name.assign(item,(unsigned int)(s_end_left - item + 1)).back() = 0;
                cimg::strellipsize(name,80,true);
                if (*title) {
                  cimg::strellipsize(title,80,true);
                  error(true,images,0,0,
                        "Operator '%s=' on variables '%s': Invalid variable name '%s'.",
                        s_operation,name.data(),title);
                } else
                  error(true,images,0,0,
                        "Operator '%s=' on variables '%s': Invalid sequence of variable names.",
                        s_operation,name.data());
              }

              // Split list of variable values in 'varvalues'.
              s = s_equal + 1;
              if (is_cond || !*s)
                varvalues.insert(1).back().assign(s,(unsigned int)std::strlen(s) + 1,1,1,1,true);
              else {
                const char *const s_end = item + std::strlen(item);
                while (s<s_end) {
                  const char *ns = std::strchr(s,',');
                  if (!ns) ns = s_end;
                  CImg<char>(s,(unsigned int)(ns - s + 1)).move_to(name);
                  name.back() = 0;
                  name.move_to(varvalues);
                  s = ns + 1;
                }
                if (*(s_end - 1)==',') CImg<char>(1,1,1,1,0).move_to(varvalues);
              }

              if (varvalues.width()!=1 && varvalues.width()!=varnames.width()) {
                name.assign(item,(unsigned int)(s_end_left - item + 1)).back() = 0;
                cimg::strellipsize(name,80,true);
                cimg::strellipsize(s_equal + 1,gmic_use_argx,80,true);
                error(true,images,0,0,
                      "Operator '%s=' on variable%s '%s': Right-hand side '%s' defines %s%d values for %s%d variables.",
                      s_operation,varnames.size()!=1?"s":"",name.data(),argx,
                      varvalues.width()<varnames.width()?"only ":"",varvalues.width(),
                      varvalues.width()>varnames.width()?"only ":"",varnames.width());
              }
            }

            // Assign or update values of variables.
            if (is_valid_expr) {
              const char *new_value = 0;
              if (!is_cond) { // Non-arithmetic assignment: A[,B,C]=a[,b,c]

                cimglist_for(varnames,k) {
                  const int l = k%varvalues.width();
                  new_value = set_variable(varnames[k],sep0,varvalues[l],0,variables_sizes);
                  if (is_verbose) {
                    cimg::strellipsize(varnames[k],gmic_use_argx,80,true);
                    cimg::strellipsize(new_value,gmic_use_argy,80,true);
                    const char *const s_sep = k==varnames.width() - 2?" and":",";
                    gmic_use_message;
                    cimg_snprintf(message,_message.width(),"'%s=%s'%s ",
                                  argx,argy,s_sep);
                    CImg<char>::string(message,false).move_to(varnames[k].assign());
                  }
                }
                if (is_verbose) {
                  (varnames>'x').move_to(name);
                  name[name.width() - 2] = 0;
                  cimg::strellipsize(name,80,true);
                  print(images,0,"Set %svariable%s %s.",
                        varnames.width()!=1?"":*item=='_'?"global ":"local ",
                        varnames.width()==1?"":"s",
                        name.data());
                }

              } else { // Arithmetic assignment: A[,B,C]*=a[,b,c]
                CImg<double> varvalues_d;
                const bool is_rounded = s_equal[1]=='_';
                s = s_equal + 1 + (is_rounded?1:0);

                if (cimg_sscanf(s,"%lf%c",&value,&end)==1)
                  varvalues_d.assign(1)[0] = is_rounded?gmic_round(value):value;
                else { // Evaluate right-hand side as a math expression.
                  CImg<T> &img = images.size()?images.back():CImg<T>::empty();
                  double fast_res;
                  if (img.__eval(s,fast_res)) // Try to get fast approximation for a single scalar first
                    varvalues_d.assign(1)[0] = fast_res;
                  else {
                    if (*s!='[') {
                      name.assign((unsigned int)std::strlen(s) + 4);
                      *name = '['; name[1] = ';'; name[name._width - 2] = ']'; name.back() = 0;
                      std::memcpy(name.data() + 2,s,name.width() - 4);
                    } else CImg<char>::string(s).move_to(name);
                    strreplace_fw(name);

                    try { img.eval(varvalues_d,name,0,0,0,0,&images); }
                    catch (CImgException &e) {
                      name.assign(item,(unsigned int)(s_end_left - item + 1)).back() = 0;
                      cimg::strellipsize(name,80,true);
                      const char *const e_ptr = std::strstr(e.what(),": ");
                      error(true,images,0,0,
                            "Operator '%s=' on variable%s '%s': Invalid right-hand side; %s",
                            s_operation,varnames.width()==1?"":"s",name.data(),e_ptr?e_ptr + 2:e.what());
                    }
                    if (varnames.width()>1 && varvalues_d.height()!=1 && varvalues_d.height()!=varnames.width()) {
                      name.assign(item,(unsigned int)(s_end_left - item + 1)).back() = 0;
                      cimg::strellipsize(name,80,true);
                      cimg::strellipsize(s_equal + 1,gmic_use_argx,80,true);
                      error(true,images,0,0,
                            "Operator '%s=' on variable%s '%s': "
                            "Right-hand side '%s' defines %s%d values for %s%d variables.",
                            s_operation,varnames.size()!=1?"s":"",name.data(),argx,
                            varvalues_d.height()<varnames.width()?"only ":"",varvalues_d.height(),
                            varvalues_d.height()>varnames.width()?"only ":"",varnames.width());
                    }
                  }
                }
                if (sep0==':' && varnames.width()==1) {
                  new_value = set_variable(varnames[0],'=',
                                           varvalues_d.value_string(',',0,is_rounded?"%g":"%.17g"),
                                           0,variables_sizes);
                  if (is_verbose) {
                    cimg::strellipsize(varnames[0],gmic_use_argx,80,true);
                    cimg::strellipsize(varvalues[0],gmic_use_argy,80,true);
                    cimg::strellipsize(new_value,gmic_use_argz,80,true);
                    print(images,0,"Set %svariable '%s:=%s'->'%s'.",
                          *item=='_'?"global ":"local ",
                          argx,argy,argz);
                  }
                } else {
                  if (is_verbose) cimg::strellipsize(varvalues[0],gmic_use_argy,80,true);
                  cimglist_for(varnames,k) {
                    const int l = k%varvalues_d.height();
                    const double vvd = is_rounded?gmic_round(varvalues_d[l]):varvalues_d[l];
                    new_value = set_variable(varnames[k],sep0==':'?'=':sep0,0,vvd,variables_sizes);
                    if (is_verbose) {
                      cimg::strellipsize(varnames[k],gmic_use_argx,80,true);
                      cimg::strellipsize(new_value,gmic_use_argz,80,true);
                      const char *const s_sep = k==varnames.width() - 2?" and":",";
                      gmic_use_message;
                      cimg_snprintf(message,_message.width(),"'%s%s=%s'->'%s'%s ",
                                    argx,s_operation,argy,argz,s_sep);
                      CImg<char>::string(message,false).move_to(varnames[k].assign());
                    }
                  }
                  if (is_verbose) {
                    (varnames>'x').move_to(name);
                    name[name.width() - 2] = 0;
                    cimg::strellipsize(name,80,true);
                    print(images,0,"Set %svariable%s %s.",
                          varnames.width()!=1?"":*item=='_'?"global ":"local ",
                          varnames.width()==1?"":"s",
                          name.data());
                  }
                }
              }
              continue;
            }
          }
        }
      }

      // Input.
      if (is_command_input) ++position;
      else {
        std::strcpy(command,"input");
        argument = item - (is_hyphen || is_plus?1:0);
        is_subst_arg = is_subst_item;
        *s_selection = 0;
      }
      gmic_substitute_args(true);
      if (!is_selection || !selection) selection.assign(1,1,1,1,images.size());

      CImg<char> indicesy(256), indicesz(256), indicesc(256);
      double dx = 0, dy = 1, dz = 1, dc = 1;
      int nb = 1;
      CImg<unsigned int> indx, indy, indz, indc;
      sepx = sepy = sepz = sepc = *indices = *indicesy = *indicesz = *indicesc = *argx = *argy = *argz = *argc = 0;

      CImg<char> arg_input(argument,(unsigned int)std::strlen(argument) + 1);
      strreplace_fw(arg_input);

      CImg<char> _gmic_selection;
      if (is_verbose) selection2string(selection,images_names,0,_gmic_selection);

      char *last_x = std::strrchr(arg_input,'x');
      if (last_x && cimg_sscanf(last_x + 1,"%d%c",&nb,&end)==1 && nb>0) *last_x = 0; else { last_x = 0; nb = 1; }
      unsigned int larg = 0;

      if (*arg_input=='0' && !arg_input[1]) {

        // Empty image(s).
        if (nb==1)
          print(images,0,"Input empty image at position%s",
                _gmic_selection.data());
        else
          print(images,0,"Input %u empty images at position%s",
                nb,_gmic_selection.data());
        g_list.assign(nb);
        CImg<char>::string("[empty]").move_to(g_list_c);
        if (--nb) {
          g_list_c.insert(nb);
          for (unsigned int i = 1; i<g_list_c.size(); ++i)
            g_list_c[i - 1].get_copymark().move_to(g_list_c[i]);
        }

      } else if (*arg_input=='(' && arg_input[(larg = (unsigned int)std::strlen(arg_input)) - 1]==')') {
        CImg<T> img;
        char delimiter = 0;

        if (arg_input[1]=='\'' &&
            ((larg>3 && arg_input[larg - 2]=='\'') ||
             (larg>5 && arg_input[larg - 3]==':' && arg_input[larg - 4]=='\'' &&
              ((delimiter = arg_input[larg-2])==',' || delimiter==';' || delimiter=='/' || delimiter=='^' ||
               is_xyzc(delimiter))))) {

          // String encoded as an image.
          CImg<char> str(arg_input.data() + 2,larg - (delimiter?5:3));
          str.back() = 0;
          cimg::strunescape(str);
          img.assign((unsigned char*)str.data(),(unsigned int)std::strlen(str));
          if (delimiter && delimiter!=',' && delimiter!='x')
            img.unroll(delimiter==';' || delimiter=='y'?'y':
                       delimiter=='/' || delimiter=='z'?'z':'c');
        } else {

          // New IxJxKxL image specified as array.
          unsigned int l, cx = 0, cy = 0, cz = 0, cc = 0, maxcx = 0, maxcy = 0, maxcz = 0;
          const char *nargument = 0;
          CImg<char> s_value(256);
          char separator = 0, unroll_axis = 0, permute_axes[5] = { 0 }, c;

          for (nargument = arg_input.data() + 1; *nargument; ) {
            *s_value = separator = 0;
            char *pd = s_value;
            // Do something faster than 'scanf("%255[0-9.eEinfa+-]")'.
            for (l = 0; l<255 && (((c=*nargument)>='0' && c<='9') || c=='.' || c=='e' || c=='E' || c=='i' || c=='n'
                                  || c=='f' || c=='a' || c=='+' || c=='-'); ++l) *(pd++) = *(nargument++);
            if (l<255) *pd = 0; else arg_error("input");
            if (*nargument) separator = *(nargument++);
            if ((separator=='^' || separator=='/' || separator==';' || separator==',' ||
                 separator==')' || separator==':') &&
                cimg_sscanf(s_value,"%lf%c",&value,&end)==1) {
              if (cx>maxcx) maxcx = cx;
              if (cy>maxcy) maxcy = cy;
              if (cz>maxcz) maxcz = cz;
              if (cx>=img._width || cy>=img._height || cz>=img._depth || cc>=img._spectrum)
                img.resize(cx>=img._width?7*cx/4 + 1:std::max(1U,img._width),
                           cy>=img._height?7*cy/4 + 1:std::max(1U,img._height),
                           cz>=img._depth?7*cz/4 + 1:std::max(1U,img._depth),
                           cc>=img._spectrum?7*cc/4 + 1:std::max(1U,img._spectrum),0);
              img(cx,cy,cz,cc) = (T)value;
              switch (separator) {
              case '^' : cx = cy = cz = 0; ++cc; break;
              case '/' : cx = cy = 0; ++cz; break;
              case ';' : cx = 0; ++cy; break;
              case ',' : ++cx; break;
              case ':' : {
                c = *nargument;
                if ((is_xyzc(c) || c==',' || c==';' || c=='/' || c=='^') &&
                    nargument[1]==')' && !nargument[2]) { unroll_axis = c; nargument+=2; }
                else if (is_xyzc(c) && is_xyzc(nargument[1]) && is_xyzc(nargument[2]) && is_xyzc(nargument[3]) &&
                         nargument[4]==')' && !nargument[5]) {
                  std::memcpy(permute_axes,nargument,4);
                  nargument+=5;
                }
                else arg_error("input");
              } break;
              case ')' : break;
              default : arg_error("input");
              }
            } else arg_error("input");
          }
          img.resize(maxcx + 1,maxcy + 1,maxcz + 1,cc + 1,0);
          if (unroll_axis) img.unroll(unroll_axis=='x' || unroll_axis==','?'x':
                                      unroll_axis=='y' || unroll_axis==';'?'y':
                                      unroll_axis=='z' || unroll_axis=='/'?'z':'c');
          else if (*permute_axes) img.permute_axes(permute_axes);
          print(images,0,"Input image at position%s, with values %s",
                _gmic_selection.data(),
                gmic_argument_text_printed());
        }
        img.move_to(g_list);
        CImg<char>::string(arg_input).move_to(g_list_c);
        if (--nb) {
          g_list.insert(nb,g_list[0]);
          g_list_c.insert(nb);
          for (unsigned int i = 1; i<g_list_c.size(); ++i)
            g_list_c[i - 1].get_copymark().move_to(g_list_c[i]);
        }

      } else if (*arg_input==gmic_store &&
                 cimg_sscanf(arg_input.data() + 1,"*store/%255[a-zA-Z0-9_]%c",&(*gmic_use_argx=0),&end)==1 &&
                 (*argx<'0' || *argx>'9')) {
        if (last_x) *last_x = 'x'; // Restore full input argument

        // Image-encoded variable.
        print(images,0,
              "Input image from variable '%s', at position%s",
              argx,_gmic_selection.data());
        unsigned int l_value = 0;
        const CImg<char> svalue = get_variable(argx,0,0,&l_value);
        try {
          if (!svalue) throw CImgArgumentException(0);
          CImgList<T>::get_unserialize(svalue,l_value + 1).move_to(g_list);
        } catch (CImgArgumentException&) {
          error(true,images,0,0,
                "Command 'input': Variable '%s' has not been assigned with command 'store'.",
                argx);
        }
        if (g_list.width()==1 && !g_list_c) // Empty list
          g_list.assign();
        else { // Non-empty list
          g_list_c.assign();
          const CImg<T> &arg = g_list.back();
          const unsigned int pend = (unsigned int)arg.size();
          for (unsigned int p = 4; p<pend; ) { // Retrieve list of image names
            unsigned int np = p;
            while (np<pend && arg[np]) ++np;
            if (np<pend) CImg<T>(arg.data(p),1,++np - p,1,1,true).move_to(g_list_c);
            p = np;
          }
          cimglist_for(g_list_c,q) g_list_c[q].unroll('x');
          if (g_list_c.size()!=g_list.size() - 1)
            error(true,images,0,0,
                  "Command 'input': Invalid binary encoding of variable '%s' "
                  "(%d items, %d names)",
                  argx,(int)g_list.size() - 1,(int)g_list_c.size());
          g_list.remove();
        }
      } else if (cimg_sscanf(arg_input,"[%255[a-zA-Z_0-9%.eE%^,:+-]%c%c",gmic_use_indices,&sep,&end)==2 && sep==']') {

        // Nb copies of existing image(s).
        const CImg<unsigned int> inds = selection2cimg(indices,images.size(),images_names,"input");
        CImg<char> s_tmp;
        if (is_verbose) selection2string(inds,images_names,1,s_tmp);
        if (nb!=1)
          print(images,0,"Input %u copies of image%s at position%s",
                (unsigned int)nb,
                s_tmp.data(),
                _gmic_selection.data());
        else
          print(images,0,"Input copy of image%s at position%s",
                s_tmp.data(),
                _gmic_selection.data());

        for (int i = 0; i<nb; ++i) cimg_foroff(inds,l) {
            g_list.insert(gmic_check(images[inds[l]]));
            (i?g_list_c[l + (i - 1)*inds.height()]:images_names[inds[l]]).get_copymark().move_to(g_list_c);
          }

      } else if ((sep=0,true) &&
                 (cimg_sscanf(arg_input,"%255[][a-zA-Z0-9_.eE%+-]%c",
                              gmic_use_argx,&end)==1 ||
                  cimg_sscanf(arg_input,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-]%c",
                              argx,gmic_use_argy,&end)==2 ||
                  cimg_sscanf(arg_input,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],"
                              "%255[][a-zA-Z0-9_.eE%+-]%c",
                              argx,argy,gmic_use_argz,&end)==3 ||
                  cimg_sscanf(arg_input,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],"
                              "%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-]%c",
                              argx,argy,argz,gmic_use_argc,&end)==4 ||
                  cimg_sscanf(arg_input,"%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],"
                              "%255[][a-zA-Z0-9_.eE%+-],%255[][a-zA-Z0-9_.eE%+-],%c",
                              argx,argy,argz,argc,&sep)==5) &&
                 ((cimg_sscanf(argx,"[%255[a-zA-Z0-9_.%+-]%c%c",gmic_use_indices,&sepx,&end)==2 &&
                   sepx==']' &&
                   (indx=selection2cimg(indices,images.size(),images_names,"input")).height()==1) ||
                  (cimg_sscanf(argx,"%lf%c",&dx,&end)==1 && dx>=1) ||
                  (cimg_sscanf(argx,"%lf%c%c",&dx,&sepx,&end)==2 && dx>0 && sepx=='%')) &&
                 (!*argy ||
                  (cimg_sscanf(argy,"[%255[a-zA-Z0-9_.%+-]%c%c",indicesy.data(),&sepy,&end)==2 &&
                   sepy==']' &&
                   (indy=selection2cimg(indicesy,images.size(),images_names,"input")).height()==1) ||
                  (cimg_sscanf(argy,"%lf%c",&dy,&end)==1 && dy>=1) ||
                  (cimg_sscanf(argy,"%lf%c%c",&dy,&sepy,&end)==2 && dy>0 && sepy=='%')) &&
                 (!*argz ||
                  (cimg_sscanf(argz,"[%255[a-zA-Z0-9_.%+-]%c%c",indicesz.data(),&sepz,&end)==2 &&
                   sepz==']' &&
                   (indz=selection2cimg(indicesz,images.size(),images_names,"input")).height()==1) ||
                  (cimg_sscanf(argz,"%lf%c",&dz,&end)==1 && dz>=1) ||
                  (cimg_sscanf(argz,"%lf%c%c",&dz,&sepz,&end)==2 && dz>0 && sepz=='%')) &&
                 (!*argc ||
                  (cimg_sscanf(argc,"[%255[a-zA-Z0-9_.%+-]%c%c",indicesc.data(),&sepc,&end)==2 &&
                   sepc==']' &&
                   (indc=selection2cimg(indicesc,images.size(),images_names,"input")).height()==1) ||
                  (cimg_sscanf(argc,"%lf%c",&dc,&end)==1 && dc>=1) ||
                  (cimg_sscanf(argc,"%lf%c%c",&dc,&sepc,&end)==2 && dc>0 && sepc=='%'))) {

        // New image with specified dimensions and optionally values.
        if (indx) { dx = (float)gmic_check(images[*indx]).width(); sepx = 0; }
        if (indy) { dy = (float)gmic_check(images[*indy]).height(); sepy = 0; }
        if (indz) { dz = (float)gmic_check(images[*indz]).depth(); sepz = 0; }
        if (indc) { dc = (float)gmic_check(images[*indc]).spectrum(); sepc = 0; }
        int idx = 0, idy = 0, idz = 0, idc = 0;
        const CImg<T>& img = images.size()?gmic_check(images.back()):CImg<T>::empty();
        if (sepx=='%') { idx = (int)cimg::round(dx*img.width()/100); if (!idx) ++idx; }
        else idx = (int)cimg::round(dx);
        if (sepy=='%') { idy = (int)cimg::round(dy*img.height()/100); if (!idy) ++idy; }
        else idy = (int)cimg::round(dy);
        if (sepz=='%') { idz = (int)cimg::round(dz*img.depth()/100); if (!idz) ++idz; }
        else idz = (int)cimg::round(dz);
        if (sepc=='%') { idc = (int)cimg::round(dc*img.spectrum()/100); if (!idc) ++idc; }
        else idc = (int)cimg::round(dc);
        if (idx<=0 || idy<=0 || idz<=0 || idc<=0) arg_error("input");
        CImg<char> s_values;
        if (sep) {
          const char *_s_values = arg_input.data() + std::strlen(argx) + std::strlen(argy) +
            std::strlen(argz) + std::strlen(argc) + 4;
          s_values.assign(_s_values,(unsigned int)std::strlen(_s_values) + 1);
          cimg::strpare(s_values,'\'',true,false);
          strreplace_fw(s_values);
          CImg<char> s_values_text(72);
          const unsigned int l = (unsigned int)std::strlen(s_values);
          if (l>=72) {
            std::memcpy(s_values_text.data(),s_values.data(),32);
            std::memcpy(s_values_text.data() + 32,"(...)",5);
            std::memcpy(s_values_text.data() + 37,s_values.data() + l - 34,35);  // Last '\0' is included
          } else std::strcpy(s_values_text,s_values);

          if (nb==1)
            print(images,0,"Input image at position%s, with values '%s'",
                  _gmic_selection.data(),s_values_text.data());
          else
            print(images,0,"Input %u images at position%s, with values '%s'",
                  nb,_gmic_selection.data(),s_values_text.data());
        } else
          print(images,0,"Input black image at position%s",
                _gmic_selection.data());
        CImg<T> new_image(idx,idy,idz,idc);
        if (s_values) {
          new_image.fill(s_values.data(),true,true,&images);
          gmic_use_title;
          cimg_snprintf(title,_title.width(),"[%s]",s_values.data());
          CImg<char>::string(title).move_to(g_list_c);
        } else { new_image.fill((T)0); CImg<char>::string("[unnamed]").move_to(g_list_c); }
        new_image.move_to(g_list);
        if (--nb) {
          g_list.insert(nb,g_list[0]);
          g_list_c.insert(nb);
          for (unsigned int i = 1; i<g_list_c.size(); ++i)
            g_list_c[i - 1].get_copymark().move_to(g_list_c[i]);
        }

      } else {

        // Input filename.
        char cext[12];
        CImg<char> _filename(4096), filename_tmp(256), options(256);
        *cext = *_filename = *filename_tmp = *options = 0;
        bool is_network_file = false;
        if (cimg_sscanf(argument,"%11[a-zA-Z0-9]:%4095[^,],%255s",
                        cext,_filename.data(),options.data())<2 ||
            !cext[1] || // Length of 'ext' must be >=2 (avoid case 'C:\\...' on Windows)
            !cimg::strcasecmp(cext,"http") || !cimg::strcasecmp(cext,"https")) {
          *cext = *_filename = *options = 0;
          if (cimg_sscanf(argument,"%4095[^,],%255s",_filename.data(),options.data())!=2) {
            std::strncpy(_filename,argument,_filename.width() - 1);
            _filename[_filename.width() - 1] = 0;
          }
        }
        strreplace_fw(_filename);
        strreplace_fw(options);
        CImg<char> _filename0 = CImg<char>::string(_filename);
        const char *const filename0 = _filename0.data();

        // Test for network file requests.
        if (!cimg::strncasecmp(_filename,"http://",7) ||
            !cimg::strncasecmp(_filename,"https://",8)) {
          try {
            cimg::load_network(_filename,filename_tmp,network_timeout,true,0,"gmic");
          } catch (CImgIOException&) {
            print(images,0,"Input file '%s' at position%s",
                  filename0,
                  _gmic_selection.data());
            error(true,images,0,0,
                  "Unreachable network file '%s'.",
                  gmic_argument_text());
          }
          is_network_file = true;
          std::strncpy(_filename,filename_tmp,_filename.width() - 1);
          _filename[_filename.width() - 1] = 0;
          *filename_tmp = 0;
        }

        if (*cext) { // Force input to be read as a '.ext' file : generate random filename
          if (*_filename=='-' && (!_filename[1] || _filename[1]=='.')) {
            // Simplify filename 'ext:-.foo' as '-.ext'.
            cimg_snprintf(_filename,_filename.width(),"-.%s",cext);
            *cext = 0;
          } else {
            std::FILE *file = 0;
            do {
              cimg_snprintf(filename_tmp,filename_tmp.width(),"%s%c%s.%s",
                            cimg::temporary_path(),cimg_file_separator,
                            cimg::filenamerand(),cext);
              if ((file=cimg::std_fopen(filename_tmp,"rb"))!=0) cimg::fclose(file);
            } while (file);

            // Make a temporary copy (or link) of the original file.
#if cimg_OS==1
            const char *const _filename_path = realpath(_filename,0);
            if (!_filename_path || symlink(_filename_path,filename_tmp))
              CImg<unsigned char>::get_load_raw(_filename).save_raw(filename_tmp);
            if (_filename_path) std::free((void*)_filename_path);
#else // #if cimg_OS==1
            CImg<unsigned char>::get_load_raw(_filename).save_raw(filename_tmp);
#endif // #if cimg_OS==1
          }
        }

        const char
          *const filename = *filename_tmp?filename_tmp:_filename,
          *const ext = cimg::split_filename(filename);
        const bool is_stdin = *filename=='-' && (!filename[1] || filename[1]=='.');
        CImg<char> uext = CImg<char>::string(ext);
        cimg::lowercase(uext);

        const char *file_type = 0;
        std::FILE *const file = is_stdin?0:cimg::std_fopen(filename,"rb");
        longT _siz = 0;
        if (file) {
          std::fseek(file,0,SEEK_END);
          _siz = std::ftell(file);
          std::rewind(file);
          file_type = *uext?0:cimg::ftype(file,0);
          cimg::fclose(file);
        }
        if (!is_stdin && file && _siz==0) { // Empty file -> Insert an empty image
          g_list.assign(1);
          _filename0.move_to(g_list_c);
        } else if (!std::strcmp("off",uext) || (file_type && !std::strcmp(file_type,"off"))) {

          // 3D object .off file.
          print(images,0,"Input 3D object '%s' at position%s",
                filename0,_gmic_selection.data());

          if (*options)
            error(true,images,0,0,
                  "Command 'input': File '%s', format does not take any input options (options '%s' specified).",
                  filename0,options.data());

          CImg<float>::get_load_off(primitives,g_list_f,filename).move_to(vertices);
          const CImg<float> opacities(1,primitives.size(),1,1,1);
          vertices.object3dtoCImg3d(primitives,g_list_f,opacities,false).move_to(g_list);
          primitives.assign();
          g_list_f.assign();
          _filename0.move_to(g_list_c);
        } else if (!std::strcmp(uext,"cimg") && *options) {

          // Part of a .cimg file (non-compressed).
          double
            n0 = -1, x0 = -1, y0 = -1, z0 = -1, c0 = -1,
            n1 = -1, x1 = -1, y1 = -1, z1 = -1, c1 = -1;
          if ((cimg_sscanf(options,"%lf,%lf%c",
                           &n0,&n1,&end)==2 ||
               cimg_sscanf(options,"%lf,%lf,%lf,%lf%c",
                           &n0,&n1,&x0,&x1,&end)==4 ||
               cimg_sscanf(options,"%lf,%lf,%lf,%lf,%lf,%lf%c",
                           &n0,&n1,&x0,&y0,&x1,&y1,&end)==6 ||
               cimg_sscanf(options,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c",
                           &n0,&n1,&x0,&y0,&z0,&x1,&y1,&z1,&end)==8 ||
               cimg_sscanf(options,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c",
                           &n0,&n1,&x0,&y0,&z0,&c0,&x1,&y1,&z1,&c1,&end)==10) &&
              (n0==-1 || n0>=0) && (n1==-1 || n1>=0) &&
              (x0==-1 || x0>=0) && (x1==-1 || x1>=0) &&
              (y0==-1 || y0>=0) && (y1==-1 || y1>=0) &&
              (z0==-1 || z0>=0) && (z1==-1 || z1>=0) &&
              (c0==-1 || c0>=0) && (c1==-1 || c1>=0)) {
            n0 = cimg::round(n0); n1 = cimg::round(n1); x0 = cimg::round(x0); x1 = cimg::round(x1);
            y0 = cimg::round(y0); y1 = cimg::round(y1); z0 = cimg::round(z0); z1 = cimg::round(z1);
            c0 = cimg::round(c0); c1 = cimg::round(c1);
            if (c0==-1 && c1==-1) {
              if (z0==-1 && z1==-1) {
                if (y0==-1 && y1==-1) {
                  if (x0==-1 && x1==-1) {
                    print(images,0,"Input crop [%d] -> [%d] of file '%s' at position%s",
                          (int)n0,(int)n1,
                          filename0,_gmic_selection.data());
                    g_list.load_cimg(filename,
                                     (unsigned int)n0,(unsigned int)n1,
                                     0U,0U,0U,0U,~0U,~0U,~0U,~0U);
                  } else {
                    print(images,0,"Input crop [%d](%d) -> [%d](%d) of file '%s' at position%s",
                          (int)n0,(int)x0,(int)n1,(int)x1,
                          filename0,_gmic_selection.data());
                    g_list.load_cimg(filename,
                                     (unsigned int)n0,(unsigned int)n1,
                                     (unsigned int)x0,0U,0U,0U,
                                     (unsigned int)x1,~0U,~0U,~0U);
                  }
                } else {
                  print(images,0,"Input crop [%d](%d,%d) -> [%d](%d,%d) of file '%s' at position%s",
                        (int)n0,(int)n1,(int)x0,(int)y0,(int)x1,(int)y1,
                        filename0,_gmic_selection.data());
                  g_list.load_cimg(filename,
                                   (unsigned int)n0,(unsigned int)n1,
                                   (unsigned int)x0,(unsigned int)y0,0U,0U,
                                   (unsigned int)x1,(unsigned int)y1,~0U,~0U);
                }
              } else {
                print(images,0,"Input crop [%d](%d,%d,%d) -> [%d](%d,%d,%d) of file '%s' "
                      "at position%s",
                      (int)n0,(int)n1,(int)x0,(int)y0,(int)z0,(int)x1,(int)y1,(int)z1,
                      filename0,_gmic_selection.data());
                g_list.load_cimg(filename,
                                 (unsigned int)n0,(unsigned int)n1,
                                 (unsigned int)x0,(unsigned int)y0,(unsigned int)z0,0U,
                                 (unsigned int)x1,(unsigned int)y1,(unsigned int)z1,~0U);
              }
            } else {
              print(images,0,"Input crop [%d](%d,%d,%d,%d) -> [%d](%d,%d,%d,%d) of file '%s' "
                    "at position%s",
                    (int)n0,(int)n1,
                    (int)x0,(int)y0,(int)z0,(int)c0,
                    (int)x1,(int)y1,(int)z1,(int)c1,
                    filename0,_gmic_selection.data());
              g_list.load_cimg(filename,
                               (unsigned int)n0,(unsigned int)n1,
                               (unsigned int)x0,(unsigned int)y0,
                               (unsigned int)z0,(unsigned int)c0,
                               (unsigned int)x1,(unsigned int)y1,
                               (unsigned int)z1,(unsigned int)c1);
            }

            if (g_list) {
              _filename0.move_to(g_list_c);
              if (g_list.size()>1) {
                g_list_c.insert(g_list.size() - 1);
                for (unsigned int i = 1; i<g_list_c.size(); ++i)
                  g_list_c[i - 1].get_copymark().move_to(g_list_c[i]);
              }
            }
          } else
            error(true,images,0,0,
                  "Command 'input': .cimg file '%s', invalid file options '%s'.",
                  filename0,options.data());

        } else if (!std::strcmp(uext,"gmz")) {
          print(images,0,"Input file '%s' at position%s",
                filename0,
                _gmic_selection.data());
          g_list.load_cimg(filename);
          bool is_gmz = false;
          const CImg<char> back = g_list?CImg<char>(g_list.back()):CImg<char>::empty();
          if (back.width()==1 && back.depth()==1 && back.spectrum()==1 &&
              back[0]=='G' && back[1]=='M' && back[2]=='Z' && !back[3]) {
            const unsigned int pend = (unsigned int)back.size();
            for (unsigned int p = 4; p<pend; ) { // Retrieve list of image names
              unsigned int np = p;
              while (np<back._height && back[np]) ++np;
              if (np<back._height) CImg<char>(back.data(p),1,++np - p,1,1,true).move_to(g_list_c);
              p = np;
            }
            if (g_list_c.size()==g_list.size()-1) {
              is_gmz = true;
              cimglist_for(g_list_c,l) g_list_c[l].unroll('x');
              g_list.remove();
            }
          }
          if (!is_gmz)
            error(true,images,0,0,"Command 'input': File '%s' is not in .gmz format (magic number not found).",
                  filename0);
          if (g_list.size()!=g_list_c.size())
            error(true,images,0,0,"Command 'input': File '%s' is not in .gmz format "
                  "(numbers of images and names do not match).",
                  filename0);

        } else if (!std::strcmp(uext,"cimg") || !std::strcmp(uext,"cimgz")) {
          print(images,0,"Input file '%s' at position%s",
                filename0,
                _gmic_selection.data());
          g_list.load_cimg(filename);
          if (g_list) {
            _filename0.move_to(g_list_c);
            if (g_list.size()>1) {
              g_list_c.insert(g_list.size() - 1);
              for (unsigned int i = 1; i<g_list_c.size(); ++i)
                g_list_c[i - 1].get_copymark().move_to(g_list_c[i]);
            }
          }

        } else if (!std::strcmp(uext,"avi") ||
                   !std::strcmp(uext,"mov") ||
                   !std::strcmp(uext,"asf") ||
                   !std::strcmp(uext,"divx") ||
                   !std::strcmp(uext,"flv") ||
                   !std::strcmp(uext,"mpg") ||
                   !std::strcmp(uext,"m1v") ||
                   !std::strcmp(uext,"m2v") ||
                   !std::strcmp(uext,"m4v") ||
                   !std::strcmp(uext,"mjp") ||
                   !std::strcmp(uext,"mp4") ||
                   !std::strcmp(uext,"mkv") ||
                   !std::strcmp(uext,"mpe") ||
                   !std::strcmp(uext,"movie") ||
                   !std::strcmp(uext,"ogm") ||
                   !std::strcmp(uext,"ogg") ||
                   !std::strcmp(uext,"qt") ||
                   !std::strcmp(uext,"rm") ||
                   !std::strcmp(uext,"vob") ||
                   !std::strcmp(uext,"webm") ||
                   !std::strcmp(uext,"wmv") ||
                   !std::strcmp(uext,"xvid") ||
                   !std::strcmp(uext,"mpeg")) {

          // Image sequence file.
          double first_frame = 0, last_frame = -1, step = 1;
          if ((cimg_sscanf(options,"%lf,%lf%c",&first_frame,&last_frame,&end)==2 ||
              cimg_sscanf(options,"%lf,%lf,%lf%c",&first_frame,&last_frame,&step,&end)==3) &&
              first_frame>=0 && (last_frame>=first_frame || last_frame==-1) && step>=0) {

            // Read several frames.
            step = cimg::round(step);
            const unsigned int
              _first_frame = (unsigned int)first_frame,
              _last_frame = last_frame>=0?(unsigned int)last_frame:~0U;
            if (_last_frame!=~0U)
              print(images,0,"Input frames %u...%u:%g of file '%s' at position%s",
                    _first_frame,_last_frame,step,
                    filename0,
                    _gmic_selection.data());
            else
              print(images,0,"Input frames %u...(last):%g of file '%s' at position%s",
                    _first_frame,step,
                    filename0,
                    _gmic_selection.data());

            g_list.load_video(filename,_first_frame,_last_frame,(unsigned int)step);
          } else if (cimg_sscanf(options,"%lf%c",&first_frame,&end)==1 &&
                     first_frame>=0) {
            // Read a single frame.
            const unsigned int _first_frame = (unsigned int)first_frame;
            print(images,0,"Input frame %u of file '%s' at position%s",
                  _first_frame,filename0,
                  _gmic_selection.data());
            g_list.load_video(filename,_first_frame,_first_frame);
          } else if (!*options) {
            // Read all frames.
            print(images,0,"Input all frames of file '%s' at position%s",
                  filename0,
                  _gmic_selection.data());
            g_list.load_video(filename);
          } else
            error(true,images,0,0,
                  "Command 'input': Video file '%s', invalid file options '%s'.",
                  filename0,options.data());
          if (g_list) {
            g_list_c.assign(g_list.size());
            _filename0.move_to(g_list_c[0]);
            if (g_list.size()>1)
              for (unsigned int i = 1; i<g_list.size(); ++i)
                g_list_c[i - 1].get_copymark().move_to(g_list_c[i]);
          }
        } else if (!std::strcmp(uext,"raw")) {

          // Raw file.
          dx = 0; dy = dz = dc = 1;
          cimg_uint64 offset = 0;
          *argx = 0;
          if (!*options ||
              cimg_sscanf(options,"%lf%c",&dx,&end)==1 ||
              cimg_sscanf(options,"%lf,%lf%c",&dx,&dy,&end)==2 ||
              cimg_sscanf(options,"%lf,%lf,%lf%c",&dx,&dy,&dz,&end)==3 ||
              cimg_sscanf(options,"%lf,%lf,%lf,%lf%c",&dx,&dy,&dz,&dc,&end)==4 ||
              cimg_sscanf(options,"%lf,%lf,%lf,%lf," cimg_fuint64 "%c",&dx,&dy,&dz,&dc,&offset,&end)==5 ||
              cimg_sscanf(options,"%255[a-z123468]%c",gmic_use_argx,&end)==1 ||
              cimg_sscanf(options,"%255[a-z123468],%lf%c",argx,&dx,&end)==2 ||
              cimg_sscanf(options,"%255[a-z123468],%lf,%lf%c",argx,&dx,&dy,&end)==3 ||
              cimg_sscanf(options,"%255[a-z123468],%lf,%lf,%lf%c",argx,&dx,&dy,&dz,&end)==4 ||
              cimg_sscanf(options,"%255[a-z123468],%lf,%lf,%lf,%lf%c",argx,&dx,&dy,&dz,&dc,&end)==5 ||
              cimg_sscanf(options,"%255[a-z123468],%lf,%lf,%lf,%lf," cimg_fuint64 "%c",argx,&dx,&dy,&dz,&dc,&offset,
                          &end)==6) {
            const char *const stype = *argx?argx:cimg::type<T>::string();
            dx = cimg::round(dx);
            dy = cimg::round(dy);
            dz = cimg::round(dz);
            dc = cimg::round(dc);
            if (dx<0 || dy<=0 || dz<=0 || dc<=0)
              error(true,images,0,0,
                    "Command 'input': Raw file '%s', invalid specified "
                    "dimensions %gx%gx%gx%g.",
                    filename0,dx,dy,dz,dc);

            if (offset)
              print(images,0,"Input raw file '%s' (offset: %lu) with type '%s' at position%s",
                    filename0,offset,stype,
                    _gmic_selection.data());
            else
              print(images,0,"Input raw file '%s' with type '%s' at position%s",
                    filename0,stype,
                    _gmic_selection.data());

#define gmic_load_raw(svalue_type,value_type) \
            if (!cimg::strcasecmp(stype,svalue_type)) \
              CImg<value_type>::get_load_raw(filename, \
                                             (unsigned int)dx,(unsigned int)dy, \
                                             (unsigned int)dz,(unsigned int)dc,false,false,\
                                             (cimg_ulong)offset).move_to(g_list);
            gmic_load_raw("bool",bool)
            else gmic_load_raw("uint8",cimg_uint8)
              else gmic_load_raw("int8",cimg_int8)
                else gmic_load_raw("uint16",cimg_uint16)
                  else gmic_load_raw("int16",cimg_int16)
                    else gmic_load_raw("uint32",cimg_uint32)
                      else gmic_load_raw("int32",cimg_int32)
                        else gmic_load_raw("uint64",cimg_uint64)
                          else gmic_load_raw("int64",cimg_int64)
                            else gmic_load_raw("float32",cimg_float32)
                              else gmic_load_raw("double64",cimg_float64)
                                else error(true,images,0,0,
                                           "Command 'input': Raw file '%s', "
                                           "invalid specified pixel type '%s'.\n",
                                           filename0,stype);
            _filename0.move_to(g_list_c);
          } else
            error(true,images,0,0,
                  "Command 'input': Raw file '%s', invalid file options '%s'.",
                  filename0,options.data());
        } else if (!std::strcmp(uext,"yuv")) {

          // YUV file.
          double first_frame = 0, last_frame = 0, step = 1, ch = 444;
          dx = 0; dy = 1;
          if ((err = cimg_sscanf(options,"%lf,%lf,%lf,%lf,%lf,%lf",
                                 &dx,&dy,&ch,&first_frame,&last_frame,&step))>=1) {
            dx = cimg::round(dx);
            dy = cimg::round(dy);
            const unsigned int ich = (unsigned int)cimg::round(ch);
            if (dx<=0 || dy<=0)
              error(true,images,0,0,
                    "Command 'input': YUV file '%s', specified dimensions (%g,%g) are invalid.",
                    filename0,dx,dy);
            if (ich!=420 && ich!=422 && ich!=444)
              error(true,images,0,0,
                    "Command 'input': YUV file '%s', specified chroma subsampling %g is invalid.",
                    filename0,ch);
            first_frame = cimg::round(first_frame);
            if (err>4) { // Load multiple frames
              last_frame = cimg::round(last_frame);
              step = cimg::round(step);
              print(images,0,"Input frames %g...%g:%g of YUV-%u:%u:%u file '%s' at position%s",
                    first_frame,last_frame,step,
                    ich/100,(ich/10)%10,ich%10,
                    filename0,
                    _gmic_selection.data());
              g_list.load_yuv(filename,(unsigned int)dx,(unsigned int)dy,ich,
                              (unsigned int)first_frame,(unsigned int)last_frame,
                              (unsigned int)step);
            } else if (err==4) { // Load a single frame
              print(images,0,"Input frames %g of YUV-%u:%u:%u file '%s' at position%s",
                    first_frame,
                    ich/100,(ich/10)%10,ich%10,
                    filename0,
                    _gmic_selection.data());
              g_list.load_yuv(filename,(unsigned int)dx,(unsigned int)dy,ich,
                              (unsigned int)first_frame,(unsigned int)first_frame);
            } else { // Load all frames
              print(images,0,"Input all frames of YUV-%u:%u:%u file '%s' at position%s",
                    ich/100,(ich/10)%10,ich%10,
                    filename0,
                    _gmic_selection.data());
              g_list.load_yuv(filename,(unsigned int)dx,(unsigned int)dy,(unsigned int)ch);
            }
            if (g_list) {
              _filename0.move_to(g_list_c);
              if (g_list.size()>1) {
                g_list_c.insert(g_list.size() - 1);
                for (unsigned int i = 1; i<g_list_c.size(); ++i)
                  g_list_c[i - 1].get_copymark().move_to(g_list[i]);
              }
            }
          } else
            error(true,images,0,0,
                  "Command 'input': YUV file '%s', invalid or missing file options '%s'.",
                  filename0,options.data());

        } else if (!std::strcmp(uext,"tif") || !std::strcmp(uext,"tiff") ||
                   (file_type && !std::strcmp(file_type,"tif"))) {

          // TIFF file.
          double first_frame = 0, last_frame = 0, step = 1;
          unsigned int bits_per_value = 0;
#ifdef cimg_use_tiff
          static const TIFFErrorHandler default_handler = TIFFSetWarningHandler(0);
          if (is_very_verbose) TIFFSetWarningHandler(default_handler);
          else TIFFSetWarningHandler(0);
#endif // #ifdef cimg_use_tiff
          if ((err = cimg_sscanf(options,"%lf,%lf,%lf",&first_frame,&last_frame,&step))>0) {
            first_frame = cimg::round(first_frame);
            if (err>1) { // Load multiple frames
              last_frame = cimg::round(last_frame);
              step = cimg::round(step);
              print(images,0,"Input frames %g...%g:%g of TIFF file '%s' at position%s",
                    first_frame,last_frame,step,
                    filename0,
                    _gmic_selection.data());
              g_list.load_tiff(filename,(unsigned int)first_frame,(unsigned int)last_frame,
                               (unsigned int)step,&bits_per_value);
            } else if (err==1) { // Load a single frame
              print(images,0,"Input frames %g of TIFF file '%s' at position%s",
                    first_frame,
                    filename0,
                    _gmic_selection.data());
              g_list.load_tiff(filename,(unsigned int)first_frame,(unsigned int)first_frame,1,&bits_per_value);
            }
          } else { // Load all frames
            if (*options) error(true,images,0,0,
                                "Command 'input': TIFF file '%s', "
                                "invalid file options '%s'.",
                                filename0,options.data());
            print(images,0,"Input all frames of TIFF file '%s' at position%s",
                  filename0,
                  _gmic_selection.data());
            g_list.load_tiff(filename,0,~0U,1,&bits_per_value);
          }
          if (g_list) {
            _filename0.move_to(g_list_c);
            if (g_list.size()>1) {
              g_list_c.insert(g_list.size() - 1);
              for (unsigned int i = 1; i<g_list_c.size(); ++i)
                g_list_c[i - 1].get_copymark().move_to(g_list_c[i]);
            }
          }
          cimg_snprintf(gmic_use_argx,256,"%d",bits_per_value);
          CImg<char>::string(argx).move_to(status);

        } else if (!std::strcmp(uext,"png")) {
          unsigned int bits_per_value = 0;
          print(images,0,"Input file '%s' at position%s",
                filename0,
                _gmic_selection.data());

          CImg<T>::get_load_png(filename,&bits_per_value).move_to(g_list);
          _filename0.move_to(g_list_c);
          cimg_snprintf(gmic_use_argx,256,"%d",bits_per_value);
          CImg<char>::string(argx).move_to(status);

        } else if (!std::strcmp(uext,"pdf")) {
          float resolution = 400;
          if (!*options || cimg_sscanf(options,"%f%c",&resolution,&end)==1) {
            const unsigned int _resolution = (int)cimg::round(std::max(resolution,20.0f));
            print(images,0,"Input file '%s' at position%s, with resolution %u",
                  filename0,_gmic_selection.data(),_resolution);
            _filename0.move_to(g_list_c);
            CImg<T>::get_load_pdf_external(filename,_resolution).move_to(g_list);
          } else
            error(true,images,0,0,
                  "Command 'input': File '%s', "
                  "invalid file options '%s'.",
                  filename0,options.data());

        } else if ((allow_main_ && !*uext) || !std::strcmp(uext,"gmic")) {

          // G'MIC command file.
          const bool add_debug_info = (*options!='0');
          print(images,0,"Input custom command file '%s'%s",
                filename0,!add_debug_info?" without debug info":"");
          unsigned int count_new = 0, count_replaced = 0;
          std::FILE *const gfile = cimg::fopen(filename,"rb");

          bool is_main_ = false, is_add_error = false;
          status.move_to(o_status); // Save status because 'add_commands' can change it, with error()
          int o_verbosity = verbosity;
          const bool o_is_debug = is_debug;
          verbosity = 0;
          is_debug = false;
          try {
            add_commands(gfile,filename,o_is_debug,&count_new,&count_replaced,
                         allow_main_ && callstack.size()==1 && !is_command_input?&is_main_:0);
          } catch (...) {
            is_add_error = true; is_main_ = false;
          }
          is_debug = o_is_debug;
          verbosity = o_verbosity;
          o_status.move_to(status);
          if (is_add_error) {
            if (is_network_file)
              error(true,images,0,0,
                    "Command 'input': Unable to load custom command file '%s' from network.",
                    filename0);
            else
              error(true,images,0,0,
                    "Command 'input': File '%s' is not recognized as a custom command file.",
                    filename0);
          }
          cimg::fclose(gfile);
          if (is_verbose) {
            unsigned int count_total = 0;
            for (unsigned int l = 0; l<gmic_comslots; ++l) count_total+=commands[l].size();
            cimg::mutex(29);
            if (count_new && count_replaced)
              std::fprintf(cimg::output()," (%u new, %u replaced, total: %u).",
                           count_new,count_replaced,count_total);
            else if (count_new)
              std::fprintf(cimg::output()," (%u new, total: %u).",
                           count_new,count_total);
            else
              std::fprintf(cimg::output()," (%u replaced, total: %u).",
                           count_replaced,count_total);
            std::fflush(cimg::output());
            cimg::mutex(29,0);
          }
          if (is_main_) { // Tell parser to run '_main_' in next iteration
            verbosity++;
            --position;
            run_main_ = true;
          }
          continue;

        } else { // Other file types.

          // Check if a custom command handling requested file format exists.
          gmic_use_formula;
          cimg_snprintf(formula,_formula.width(),"+input_%s",uext.data());
          hash = hashcode(formula,false);
          if (search_sorted(formula,commands_names[hash],commands_names[hash].size(),pattern)) { // Command found
            cimg_snprintf(formula,_formula.width(),"+input_%s[] \"%s\"%s%s",uext.data(),filename0,
                          *options?",":"",options.data());
            const CImgList<char> ncommands_line = commands_line_to_CImgList(formula);
            unsigned int nposition = 0;
            CImg<char>::string("").move_to(callstack); // Anonymous scope
            _run(ncommands_line,nposition,g_list,g_list_c,images,images_names,variables_sizes,0,0,0,false);
            callstack.remove();

          } else { // Not found -> Try generic image loader

            print(images,0,"Input file '%s' at position%s",
                  filename0,
                  _gmic_selection.data());
            if (*options)
              error(true,images,0,0,
                    "Command 'input': File '%s', format does not take any input options (options '%s' specified).",
                    filename0,options.data());

            try {
              try {
                g_list.load(filename);
              } catch (CImgIOException&) {
                if (is_network_file)
                  error(true,images,0,0,
                        "Command 'input': Unable to load image file '%s' from network.",
                        filename0);
                else throw;
              }

              // If .gmz file without extension, process images names anyway.
              bool is_gmz = false;
              const CImg<char> back = g_list?CImg<char>(g_list.back()):CImg<char>::empty();
              if (back.width()==1 && back.depth()==1 && back.spectrum()==1 &&
                  back[0]=='G' && back[1]=='M' && back[2]=='Z' && !back[3]) {
                g_list_c.assign();
                const unsigned int pend = (unsigned int)back.size();
                for (unsigned int p = 4; p<pend; ) { // Retrieve list of image names
                  unsigned int np = p;
                  while (np<back._height && back[np]) ++np;
                  if (np<back._height) CImg<char>(back.data(p),1,++np - p,1,1,true).move_to(g_list_c);
                  p = np;
                }
                if (g_list_c) {
                  is_gmz = true;
                  cimglist_for(g_list_c,l) g_list_c[l].unroll('x');
                  g_list.remove();
                }
              }

              if (g_list && !is_gmz) {
                g_list_c.insert(_filename0);
                if (g_list.size()>1) {
                  g_list_c.insert(g_list.size() - 1);
                  for (unsigned int i = 1; i<g_list_c.size(); ++i)
                    g_list_c[i - 1].get_copymark().move_to(g_list_c[i]);
                }
              }
            } catch (CImgException&) {
              std::FILE *efile = 0;
              if (!(efile = cimg::std_fopen(filename,"r"))) {
                if (is_command_input || std::strchr(filename,'.')!=0)
                  error(true,images,0,0,
                        "Command 'input': Unknown filename '%s'.",
                        gmic_argument_text());
                else {
                  CImg<char>::string(filename).move_to(name);
                  const unsigned int foff = (*name=='+' || *name=='-');
                  const char *misspelled = 0;
                  char *const posb = std::strchr(name,'[');
                  if (posb) *posb = 0; // Discard selection from the command name
                  int dmin = 4;
                  // Look for a built-in command.
                  for (unsigned int l = 0; l<sizeof(builtin_commands_names)/sizeof(char*); ++l) {
                    const char *const c = builtin_commands_names[l];
                    if (c) {
                      const int d = levenshtein(c,name.data() + foff);
                      if (d<dmin) { dmin = d; misspelled = builtin_commands_names[l]; }
                    }
                  }
                  // Look for a custom command.
                  for (unsigned int i = 0; i<gmic_comslots; ++i)
                    cimglist_for(commands_names[i],l) {
                      const char *const c = commands_names[i][l].data();
                      const int d = levenshtein(c + (*c=='+'?1:0),name.data() + foff);
                      if (d<dmin) { dmin = d; misspelled = commands_names[i][l].data(); }
                    }
                  if (misspelled)
                    error(true,images,0,0,
                          "Unknown command or filename '%s'; did you mean '%s'?",
                          gmic_argument_text(),misspelled + (*misspelled=='+'?1:0));
                  else error(true,images,0,0,
                             "Unknown command or filename '%s'.",
                             gmic_argument_text());
                }
              } else { std::fclose(efile); throw; }
            }
          }
        }

        if (*filename_tmp) std::remove(filename_tmp); // Clean temporary file if used
        if (is_network_file) std::remove(_filename);  // Clean temporary file if network input
      }

      if (is_verbose) {
        cimg::mutex(29);
        if (g_list) {
          const unsigned int last = g_list.size() - 1;
          if (g_list.size()==1) {
            if (g_list[0].is_CImg3d(false))
              std::fprintf(cimg::output()," (%u vertices, %u primitives).",
                           cimg::float2uint((float)g_list(0,6)),
                           cimg::float2uint((float)g_list(0,7)));
            else
              std::fprintf(cimg::output()," (1 image %dx%dx%dx%d).",
                           g_list[0].width(),g_list[0].height(),
                           g_list[0].depth(),g_list[0].spectrum());
          } else
            std::fprintf(cimg::output()," (%u images [0] = %dx%dx%dx%d, %s[%u] = %dx%dx%dx%d).",
                         g_list.size(),
                         g_list[0].width(),g_list[0].height(),
                         g_list[0].depth(),g_list[0].spectrum(),
                         last==1?"":"(...),",last,
                         g_list[last].width(),g_list[last].height(),
                         g_list[last].depth(),g_list[last].spectrum());
        } else {
          std::fprintf(cimg::output()," (no available data).");
          g_list_c.assign();
        }
        std::fflush(cimg::output());
        cimg::mutex(29,0);
      }

      for (unsigned int l = 0, lsiz = selection.height() - 1U, off = 0; l<=lsiz; ++l) {
        uind = selection[l] + off;
        off+=g_list.size();
        if (l!=lsiz) {
          images.insert(g_list,uind);
          images_names.insert(g_list_c,uind);
        } else {
          g_list.move_to(images,uind);
          g_list_c.move_to(images_names,uind);
        }
      }

      g_list.assign();
      g_list_c.assign();
      is_change = true;
    } // End main parsing loop of _run()

    // Wait for remaining threads to finish and possibly throw exceptions from threads.
    cimglist_for(gmic_threads,k) wait_threads(&gmic_threads[k],false,(T)0);
    cimglist_for(gmic_threads,k) cimg_forY(gmic_threads[k],l) {
      const gmic_exception &e = gmic_threads(k,l).exception;
      if (e._message) error(false,images,0,0,e._message);
    }

    // Post-check global environment consistency.
    if (images_names.size()!=images.size())
      error(true,images,0,0,
            "Internal error: Images (%u) and images names (%u) have different size, "
            "at return point.",
            images_names.size(),images.size());
    if (!callstack)
      error(true,images,0,0,
            "Internal error: Empty call stack at return point.");

    // Post-check call stack consistency.
    if (!is_quit && !is_return) {
      const CImg<char>& s = callstack.back();
      if (s[0]=='*' && (s[1]=='b' || s[1]=='d' || s[1]=='i' || s[1]=='r' ||
                        (s[1]=='f' && (s[4]!='e' || !is_end_local)) ||
                        (s[1]=='l' && !is_end_local))) {
        unsigned int reference_line = ~0U;
        if (cimg_sscanf(s,"*%*[a-z]#%u",&reference_line)==1)
          error(true,images,0,0,
                "A '%s' command is missing (for '%s', line #%u), before return point.",
                s[1]=='b'?"}":s[1]=='d'?"while":s[1]=='i'?"fi":"done",
                s[1]=='b'?"{":s[1]=='d'?"do":s[1]=='i'?"if":s[1]=='r'?"repeat":
                s[1]=='f'?(s[4]!='e'?"for":"foreach"):"local",
                reference_line);
        else error(true,images,0,0,
                   "A '%s' command is missing, before return point.",
                   s[1]=='d'?"while":s[1]=='i'?"fi":"done");
      }
    }
    pop_callstack(initial_callstack_size);

    // Post-check validity of shared images.
    cimglist_for(images,l) gmic_check(images[l]);

    // Display or print result.
    if (verbosity>0 && is_change && !is_quit && !is_return && callstack.size()==1 && images) {
      const CImg<char> host = get_variable("_host");
      if (host && !std::strcmp(host,"cli")) {
        if (is_display_available) {
          CImgList<unsigned int> lselection, lselection3d;
          bool is_first3d = false;
          gmic_display_window(0).assign();
          cimglist_for(images,l) {
            const bool is_3d = images[l].is_CImg3d(false);
            if (!l) is_first3d = is_3d;
            CImg<unsigned int>::vector(l).move_to(is_3d?lselection3d:lselection);
          }

          // Prepare for '_run()' used for 3D interactive viewer.
          CImgList<char> ncommands_line;
          unsigned int nposition = 0;
          if (lselection3d) {
            gmic_use_formula;
            cimg_snprintf(formula,_formula.width(),"d3d[%s]",(lselection3d>'y').value_string().data());
            commands_line_to_CImgList(formula).move_to(ncommands_line);
          }

          if (is_first3d) {
            CImg<char>::string("").move_to(callstack); // Anonymous scope
            _run(ncommands_line,nposition,images,images_names,images,images_names,variables_sizes,0,0,0,false);
            callstack.remove();
            if (lselection) display_images(images,images_names,lselection>'y',0,false);
          } else {
            if (lselection) display_images(images,images_names,lselection>'y',0,false);
            if (lselection3d) {
              CImg<char>::string("").move_to(callstack); // Anonymous scope
              _run(ncommands_line,nposition,images,images_names,images,images_names,variables_sizes,0,0,0,false);
              callstack.remove();
            }
          }
        } else {
          CImg<unsigned int> seq(1,images.width());
          cimg_forY(seq,y) seq[y] = y;
          print_images(images,images_names,seq,true);
        }
      }
      is_change = false;
    }

    if (is_debug) debug(images,"%sExit scope '%s/'.%s\n",
                        cimg::t_bold,callstack.back().data(),cimg::t_normal);

    if (callstack.size()==1) {
      if (is_quit) {
        if (verbosity>=1 || is_debug) {
          std::fputc('\n',cimg::output());
          std::fflush(cimg::output());
        }
      } else {
        print(images,0,"End G'MIC interpreter.\n");
        is_quit = true;
      }
    }
    verbosity = starting_verbosity;

  } catch (gmic_exception&) {
    // Wait for remaining threads to finish.
    cimglist_for(gmic_threads,k) wait_threads(&gmic_threads[k],true,(T)0);
    pop_callstack(initial_callstack_size);
    throw;

  } catch (CImgAbortException &) { // Special case of abort (abort from a CImg method)
    // Wait for remaining threads to finish.
    cimglist_for(gmic_threads,k) wait_threads(&gmic_threads[k],true,(T)0);
    pop_callstack(initial_callstack_size);

    // Do the same as for a cancellation point.
    const bool is_very_verbose = verbosity>1 || is_debug;
    if (is_very_verbose) print(images,0,"Abort G'MIC interpreter (caught abort signal).");
    position = commands_line.size();
    is_change = false;
    is_quit = true;

  } catch (CImgException &e) {
    // Wait for remaining threads to finish.
    cimglist_for(gmic_threads,k) wait_threads(&gmic_threads[k],true,(T)0);
    pop_callstack(initial_callstack_size);

    const char *const e_ptr = e.what() + (!std::strncmp(e.what(),"[gmic_math_parser] ",19)?19:0);
    CImg<char> error_message(e_ptr,(unsigned int)std::strlen(e_ptr) + 1);
    const char *const s_fopen = "cimg::fopen(): Failed to open file '";
    const unsigned int l_fopen = (unsigned int)std::strlen(s_fopen);
    if (!std::strncmp(error_message,s_fopen,l_fopen) &&
        !std::strcmp(error_message.end() - 18,"' with mode 'rb'.")) {
      error_message[error_message.width() - 18] = 0;
      error(true,images,0,0,"Unknown filename '%s'.",error_message.data(l_fopen));
    }
    for (char *str = std::strstr(error_message,"CImg<"); str; str = std::strstr(str,"CImg<")) {
      str[0] = 'g'; str[1] = 'm'; str[2] = 'i'; str[3] = 'c';
    }
    for (char *str = std::strstr(error_message,"CImgList<"); str; str = std::strstr(str,"CImgList<")) {
      str[0] = 'g'; str[1] = 'm'; str[2] = 'i'; str[3] = 'c';
    }
    for (char *str = std::strstr(error_message,"cimg:"); str; str = std::strstr(str,"cimg:")) {
      str[0] = 'g'; str[1] = 'm'; str[2] = 'i'; str[3] = 'c';
    }
    if (*command) {
      const char *em = error_message.data();
      if (!std::strncmp("gmic<",em,5)) {
        em = std::strstr(em,"(): ");
        if (em) em+=4; else em = error_message.data();
      }
      error(true,images,0,command,"Command '%s': %s",command,em);
    } else error(true,images,0,0,"%s",error_message.data());
  }
  if (!is_end_local) debug_line = initial_debug_line;
  else {
    if (next_debug_line!=~0U) { debug_line = next_debug_line; next_debug_line = ~0U; }
    if (next_debug_filename!=~0U) { debug_filename = next_debug_filename; next_debug_filename = ~0U; }
  }

  // Remove/modify current run from managed list of gmic runs.
  cimg::mutex(24);
  for (int k = grl.width() - (ind_run<grl._width?0:1); k>=0; --k) {
    const int _k = k>=grl.width()?ind_run:k; // First try is 'ind_run' if defined
    CImg<void*> &_gr = grl[_k];
    if (_gr && _gr[0]==this) { if (push_new_run) grl.remove(_k); else gr.swap(_gr); break; }
  }
  cimg::mutex(24,0);

  return *this;
}

// Explicitly instantiate constructors and destructor when compiling the library.
#define gmic_export_funcs(pt) \
template gmic::gmic(const char *const commands_line, const char *const custom_commands, \
                    const bool include_stdlib, float *const p_progress, bool *const p_is_abort, \
                    const pt& pixel_type); \
template gmic::gmic(const char *const commands_line, \
                    CImgList<pt>& images, CImgList<char>& images_names, \
                    const char *const custom_commands, const bool include_stdlib, \
                    float *const p_progress, bool *const p_is_abort); \
template gmic& gmic::assign(const char *const commands_line, const char *const custom_commands, \
                            const bool include_stdlib, float *const p_progress, bool *const p_is_abort, \
                            const pt& pixel_type); \
template gmic& gmic::assign(const char *const commands_line, \
                            CImgList<pt>& images, CImgList<char>& images_names, \
                            const char *const custom_commands, const bool include_stdlib, \
                            float *const p_progress, bool *const p_is_abort); \
template gmic& gmic::run(const char *const commands_line, const pt& pixel_type); \
template gmic& gmic::run(const char *const commands_line, \
                         CImgList<pt> &images, CImgList<char> &images_names); \
template CImg<pt>& CImg<pt>::assign(const unsigned int size_x, const unsigned int size_y, \
                                    const unsigned int size_z, const unsigned int size_c); \
template CImgList<pt>& CImgList<pt>::assign(const unsigned int n)

#ifdef gmic_pixel_type
gmic_export_funcs(gmic_pixel_type);
#endif
#ifdef gmic_pixel_type2
gmic_export_funcs(gmic_pixel_type2);
#endif
template CImgList<char>::~CImgList();
template CImgList<char>& CImgList<char>::assign(const unsigned int n);
template bool gmic::search_sorted(const char *const str, const CImgList<char>& list,
                                  const unsigned int length, unsigned int &out_ind);
#endif // #ifdef cimg_plugin
