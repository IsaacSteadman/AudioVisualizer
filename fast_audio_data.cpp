#include <cstdint>
#include <limits>
#include <utility>
#include <cmath>
#include <type_traits>
#include <iostream>
#ifdef __linux__
#define FAST_AUDIO_DATA_EXPORT
#else
#define FAST_AUDIO_DATA_EXPORT __declspec(dllexport)
#endif
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279
#endif

using std::swap;
using std::sin;
using std::max;
using std::min;
using std::is_same_v;
using std::numeric_limits;
using std::pow;
using std::sqrt;
using std::log2;
using std::cout;

#define TRY_BLOCK_START try {
#define TRY_BLOCK_END } catch (...) { cout << "ERROR uncaught exception at " << __FILE__ << ':' << __LINE__ << '\n'; throw; }
#define TRY_BLOCK_END_EXTRA(extra) } catch (...) { cout << "ERROR uncaught exception at " << __FILE__ << ':' << __LINE__ << ", extra " << extra << '\n'; throw; }
#define TRY_BLOCK_END_EXTRA3(extra1, extra2, extra3) } catch (...) { cout << "ERROR uncaught exception at " << __FILE__ << ':' << __LINE__ << ", extra1 " << extra1 << " extra2 " << extra2 << " extra3 " << extra3 << '\n'; throw; }

static size_t used_fft_size = 0;
static size_t audio_framesize = 0;
static size_t audio_nchannels = 0;
static double *fft_data = nullptr;
static double *channel_temp_data = nullptr;
static struct PtRange {
  double pt0, pt1;
} *lst_pt_ranges = nullptr;
static size_t used_display_width = 0;
static size_t used_avg_keep = 0;
static size_t num_avg_samples = 0;
static size_t idx_insert_sample = 0;
static double *avg_keep_data = nullptr;

template<typename T, bool swapped>
static inline T read_swapped(const uint8_t *data_ptr) {
  union {
    T data;
    uint8_t bytes[sizeof(T)];
  } helper;
  helper.data = *(T *)data_ptr;
  if constexpr (swapped && sizeof(T) > 1) {
    for (size_t i = 0; i < sizeof(T) / 2; ++i) {
      swap(helper.bytes[i], helper.bytes[sizeof(T) - i - 1]);
    }
  }
  return helper.data;
}

static char *digits = "0123456789ABCDEF";

template<typename T, bool swapped>
static inline int convert_frames_to_doubles_swapped(const uint8_t *framebuffer) {
  if (fft_data == nullptr) {
    return -1;
  } else if (framebuffer == nullptr) {
    return -2;
  } else if (audio_nchannels == 0) {
    return -3;
  } else if (audio_framesize != sizeof(T)) {
    return -4;
  }
  const size_t len = used_fft_size * 2;
  // cout << "data hex (first 108 bytes):";
  // for (size_t i = 0; i < min(len * audio_nchannels * audio_framesize, (size_t)108); ++i) {
  //   const uint8_t b = framebuffer[i];
  //   cout << ' ';
  //   cout << digits[b >> 4];
  //   cout << digits[b & 0xF];
  // }
  // cout << '\n';
  if (audio_nchannels == 2) {
    for (size_t i = 0; i < len; ++i) {
      const uint8_t *base = framebuffer + i * audio_nchannels * audio_framesize;
      double sum = (read_swapped<T, swapped>(base) + read_swapped<T, swapped>(base + audio_framesize)) / 2.0;
      if constexpr (numeric_limits<T>::is_integer && !numeric_limits<T>::is_signed) {
        sum = 2.0 * sum / numeric_limits<T>::max() - 1.0;
      } else if constexpr (numeric_limits<T>::is_integer) {
        sum /= -numeric_limits<T>::min();
      }
      // if (i % 128 == 0) cout << "sum@" << i << " = " << sum << '\n';
      fft_data[i] = sum;
    }
  } else if (audio_nchannels == 1) {
    for (size_t i = 0; i < len; ++i) {
      double sum = read_swapped<T, swapped>(framebuffer + i * audio_nchannels * audio_framesize);
      if constexpr (numeric_limits<T>::is_integer && !numeric_limits<T>::is_signed) {
        sum = 2.0 * sum / numeric_limits<T>::max() - 1.0;
      } else if constexpr (numeric_limits<T>::is_integer) {
        sum /= -numeric_limits<T>::min();
      }
      fft_data[i] = sum;
    }
  } else {
    for (size_t i = 0; i < len; ++i) {
      const uint8_t *base = framebuffer + i * audio_nchannels * audio_framesize;
      double sum = 0.0;
      for (size_t j = 0; j < audio_nchannels; ++j,base += audio_framesize) {
        sum += read_swapped<T, swapped>(base);
      }
      sum /= (double) audio_nchannels;
      if constexpr (numeric_limits<T>::is_integer && !numeric_limits<T>::is_signed) {
        sum = 2.0 * sum / numeric_limits<T>::max() - 1.0;
      } else if constexpr (numeric_limits<T>::is_integer) {
        sum /= -numeric_limits<T>::min();
      }
      fft_data[i] = sum;
    }
  }
  return 0;
}

static inline void fft_internal(double* data, size_t nn)
{
  size_t n, mmax, m, j, istep, i;
  double wtemp, wr, wpr, wpi, wi, theta;
  double tempr, tempi;

  // reverse-binary reindexing
  n = nn<<1;
  j=1;
  for (i=1; i<n; i+=2) {
    if (j>i) {
      swap(data[j-1], data[i-1]);
      swap(data[j], data[i]);
    }
    m = nn;
    while (m>=2 && j>m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  };

  // here begins the Danielson-Lanczos section
  mmax=2;
  while (n>mmax) {
    istep = mmax<<1;
    theta = -(2*M_PI/mmax);
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m=1; m < mmax; m += 2) {
      for (i=m; i <= n; i += istep) {
        j=i+mmax;
        tempr = wr*data[j-1] - wi*data[j];
        tempi = wr * data[j] + wi*data[j-1];

        data[j-1] = data[i-1] - tempr;
        data[j] = data[i] - tempi;
        data[i-1] += tempr;
        data[i] += tempi;
      }
      wtemp=wr;
      wr += wr*wpr - wi*wpi;
      wi += wi*wpr + wtemp*wpi;
    }
    mmax=istep;
  }
}

static inline double getValueFromPtRangeValues(double *values, size_t len, double pt0, double pt1) {
  size_t ipt0 = pt0;
  size_t ipt1 = pt1;
  double fpt0 = pt0 - ipt0;
  double fpt1 = pt1 - ipt1;
  if (ipt0 != ipt1) {
    double sum = values[max(ipt0, (size_t)0)] * fpt0;
    const size_t end = min(ipt1 + 1, len);
    const size_t start = ipt0 + 1;
    for (size_t i = start; i < end; ++i) {
      sum += values[i];
    }
    sum += values[min(ipt1 + 1, len - 1)] * fpt1;
    return sum / (fpt0 + fpt1 + end - start);
  } else if (ipt0 + 1 < len && ipt0 >= 0) {
    double val0 = values[ipt0];
    double val1 = values[ipt0 + 1];
    return ((2.0 - fpt0 - fpt1) * val0 + (fpt0 + fpt1) * val1) / 2.0;
  } else if (ipt0 + 1 >= len) {
    return values[len - 1];
  } else if (ipt0 < 0) {
    return values[0];
  }
}

static inline double getValueFromPtRangeValues2(double *values, size_t len, double pt0, double pt1) {
  size_t ipt0 = pt0;
  size_t ipt1 = pt1;
  double fpt0 = pt0 - ipt0;
  double fpt1 = pt1 - ipt1;
  if (ipt1 - ipt0 > 1) {
    double sum = values[max(ipt0, (size_t)0)] * fpt0;
    const size_t end = min(ipt1 + 1, len);
    const size_t start = ipt0 + 1;
    for (size_t i = start; i < end; ++i) {
      sum += values[i];
    }
    sum += values[min(ipt1 + 1, len - 1)] * fpt1;
    return sum / (fpt0 + fpt1 + end - start);
  } else if (ipt1 - ipt0 == 1 && ipt0 + 1 < len && ipt0 >= 0) {
    double val0 = values[ipt0];
    double val1 = values[ipt1];
    double val2 = values[ipt1 + 1];
    return ((1.0 - fpt0) * val0 + (fpt0 + 1.0 - fpt1) * val1 + fpt1 * val2) / 2.0;
  } else if (ipt0 + 1 < len && ipt0 >= 0) {
    double val0 = values[ipt0];
    double val1 = values[ipt0 + 1];
    return ((2.0 - fpt0 - fpt1) * val0 + (fpt0 + fpt1) * val1) / 2.0;
  } else if (ipt1 + 1 >= len) {
    return values[len - 1];
  } else if (ipt0 < 0) {
    return values[0];
  }
}

const double NormC0 = pow(20.6, 2);
const double NormC1 = pow(12194, 2);
const double NormC2 = pow(107.7, 2);
const double NormC3 = pow(737.9, 2);
const double NormC4 = pow(10, .2);
const double NormC5 = pow(NormC1, 2);

// calculate the human threshold for hearing adjustment coefficient given frequency x in hertz
static inline double calc_hth_adjust_coef(double x) {
    if (x == 0) return 1;
    double a = x * x;
    double n = pow(a + NormC0, 2) * pow(a + NormC1, 2) * (a + NormC2) * (a + NormC3);
    double d = NormC4 * NormC5 * pow(a, 4);
    return n/d;
}

static inline double get_at_ipt(double *fft_data, size_t len, ptrdiff_t ipt) {
  if (ipt < 0) {
    return fft_data[0];
  } else if (ipt >= len) {
    return fft_data[len - 1];
  } else {
    fft_data[ipt];
  }
}

static inline double get_at_pt(double *fft_data, size_t len, double pt) {
  ptrdiff_t l_pt = pt;
  if (pt == l_pt) {
    return get_at_ipt(fft_data, len, l_pt);
  }
  double a = pt - l_pt;
  ptrdiff_t u_pt = l_pt + 1;
  double l_v = get_at_ipt(fft_data, len, l_pt);
  double u_v = get_at_ipt(fft_data, len, u_pt);
  return (1 - a) * l_v + a * u_v;
}

static inline double getValueFromPtRangeValues1(double *fft_data, size_t len, double pt0, double pt1) {
  ptrdiff_t ipt0 = pt0;
  ptrdiff_t ipt1 = pt1;
  if (ipt0 == ipt1) {
    return (get_at_pt(fft_data, len, pt0) + get_at_pt(fft_data, len, pt1)) / 2.0;
  }
  double fpt0 = pt0 - ipt0;
  double fpt1 = pt1 - ipt1;
  double sum = fft_data[max(ipt0, (ptrdiff_t)0)] * fpt0;
  const size_t end = min(ipt1 + 1, (ptrdiff_t)len);
  const size_t start = ipt0 + 1;
  for (size_t i = start; i < end; ++i) {
    sum += fft_data[i];
  }
  sum += fft_data[min(ipt1 + 1, (ptrdiff_t)(len - 1))] * fpt1;
  return sum / (fpt0 + fpt1 + end - start);
}

static inline void render_sdl_fft_data(void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
  fft_internal(fft_data, used_fft_size);
  double max_val = -9e99;
  TRY_BLOCK_START
  const size_t fft_data_len = used_fft_size * 2;
  for (size_t i = 1; i < used_fft_size; ++i) {
    const double x = fft_data[i];
    const double y = fft_data[fft_data_len - i];
    const double cur = sqrt(x * x + y * y) / calc_hth_adjust_coef(sample_rate_hz * i / fft_data_len);
    max_val = max(max_val, cur);
    fft_data[i] = cur;
  }
  TRY_BLOCK_END
  TRY_BLOCK_START
  for (size_t i = 1; i < used_fft_size; ++i) {
    fft_data[i] /= max_val;
  }
  fft_data[0] = 0.0;
  TRY_BLOCK_END
  if (used_avg_keep == 0) {
    TRY_BLOCK_START
    for (size_t x = 0; x < used_display_width; ++x) {
      double val;
      TRY_BLOCK_START
      val = min(1.0, max(0.0, getValueFromPtRangeValues2(fft_data, used_fft_size, lst_pt_ranges[x].pt0, lst_pt_ranges[x].pt1)));
      TRY_BLOCK_END_EXTRA(x)
      const double inv_val = 1.0 - val;
      uint8_t r = sr * inv_val + er * val;
      uint8_t g = sg * inv_val + eg * val;
      uint8_t b = sb * inv_val + eb * val;
      uint32_t color = b | (g << 8) | (r << 16);
      size_t inv_bar_height = height * inv_val;
      size_t y = height;
      TRY_BLOCK_START
      while (y-- > inv_bar_height) {
        ((uint32_t *)((uint8_t *)pixels + pitch * y))[x] = color;
      }
      if (inv_bar_height > 0) {
        while (y-- > 0) {
          ((uint32_t *)((uint8_t *)pixels + pitch * y))[x] = 0;
        }
      }
      TRY_BLOCK_END_EXTRA3(x, y, inv_bar_height)
    }
    TRY_BLOCK_END
  } else {
    TRY_BLOCK_START
    for (size_t x = 0; x < used_display_width; ++x) {
      double sample;
      TRY_BLOCK_START
      sample = min(1.0, max(0.0, getValueFromPtRangeValues2(fft_data, used_fft_size, lst_pt_ranges[x].pt0, lst_pt_ranges[x].pt1)));
      TRY_BLOCK_END_EXTRA(x)
      double val = sample;
      double *base_avg_data = avg_keep_data + x * used_avg_keep;
      for (size_t i = 0; i < num_avg_samples; ++i) {
        val += base_avg_data[i];
      }
      base_avg_data[idx_insert_sample] = sample;
      val /= used_avg_keep + 1;
      const double inv_val = 1.0 - val;
      uint8_t r = sr * inv_val + er * val;
      uint8_t g = sg * inv_val + eg * val;
      uint8_t b = sb * inv_val + eb * val;
      uint32_t color = b | (g << 8) | (r << 16);
      size_t inv_bar_height = height * inv_val;
      size_t y = height;
      TRY_BLOCK_START
      while (y-- > inv_bar_height) {
        ((uint32_t *)((uint8_t *)pixels + pitch * y))[x] = color;
      }
      if (inv_bar_height > 0) {
        while (y-- > 0) {
          ((uint32_t *)((uint8_t *)pixels + pitch * y))[x] = 0;
        }
      }
      TRY_BLOCK_END_EXTRA3(x, y, inv_bar_height)
    }
    idx_insert_sample = (idx_insert_sample + 1) % used_avg_keep;
    num_avg_samples = min(num_avg_samples + 1, used_avg_keep);
    TRY_BLOCK_END
  }
}

extern "C" {
  FAST_AUDIO_DATA_EXPORT void fft(double *data, size_t nn) {
    fft_internal(data, nn);
  }
  FAST_AUDIO_DATA_EXPORT void fad_destroy_buffers() {
    if (fft_data) {
      delete[] fft_data;
      fft_data = nullptr;
    }
    if (lst_pt_ranges) {
      delete[] lst_pt_ranges;
      lst_pt_ranges = nullptr;
    }
    if (channel_temp_data) {
      delete[] channel_temp_data;
      channel_temp_data = nullptr;
    }
    if (avg_keep_data) {
      delete[] avg_keep_data;
      avg_keep_data = nullptr;
      num_avg_samples = 0;
      idx_insert_sample = 0;
    }
  }
  FAST_AUDIO_DATA_EXPORT void fad_build_buffers(size_t framesize, size_t nchannels, size_t fft_size, size_t display_width) {
    fad_destroy_buffers();
    fft_data = new double[2 * fft_size];
    used_fft_size = fft_size;
    lst_pt_ranges = new PtRange[display_width];
    used_display_width = display_width;
    audio_framesize = framesize;
    audio_nchannels = nchannels;
    if (nchannels > 2) {
      channel_temp_data = new double[nchannels];
    }
    if (used_avg_keep) {
      avg_keep_data = new double[used_display_width * used_avg_keep];
    }
    const size_t start = 1;
    const size_t stop = used_fft_size;
    const double scale_start = log2(start);
    const double scale_stop = log2(stop);
    const double scale_step = (scale_stop - scale_start) / used_display_width;
    for (size_t i = 0; i < used_display_width; ++i) {
      lst_pt_ranges[i].pt0 = max(0.0, exp2((i - 0.5) * scale_step) * start);
      lst_pt_ranges[i].pt1 = min(used_fft_size - 1.0, exp2((i + 0.5) * scale_step) * start);
    }
    cout << "BUFFER fft_data (address = 0x" << (void *)fft_data << " size = " << 2 * fft_size * sizeof(double) << " bytes)\n";
    cout << "BUFFER lst_pt_ranges (address = 0x" << (void *)lst_pt_ranges << " size = " << display_width * sizeof(PtRange) << " bytes)\n";
    cout << "BUFFER channel_temp_data (address = 0x" << (void *)channel_temp_data << " size = " << nchannels * sizeof(double) << " bytes)\n";
  }
  FAST_AUDIO_DATA_EXPORT void fad_set_averaging(size_t averaging) {
    if (used_avg_keep == averaging || averaging > 100) return;
    if (avg_keep_data) {
      delete[] avg_keep_data;
      avg_keep_data = nullptr;
      num_avg_samples = 0;
      idx_insert_sample = 0;
    }
    used_avg_keep = averaging;
    if (used_avg_keep) {
      avg_keep_data = new double[used_display_width * used_avg_keep];
    }
  }
  FAST_AUDIO_DATA_EXPORT size_t fad_get_averaging() {
    return used_avg_keep;
  }
  FAST_AUDIO_DATA_EXPORT size_t fad_get_audio_buffer_size() {
    return audio_nchannels * audio_framesize * 2 * used_fft_size;
  }
  // ne means native/default endian-ness
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_int32ne_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<int32_t, false>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_int16ne_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<int16_t, false>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_uint32ne_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<uint32_t, false>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_uint16ne_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<uint16_t, false>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_float64ne_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<double, false>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_float32ne_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<float, false>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  // fe means foreign/non-native/non-default endian-ness
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_int32fe_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<int32_t, true>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_int16fe_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<int16_t, true>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_uint32fe_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<uint32_t, true>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_uint16fe_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<uint16_t, true>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_float64fe_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<double, true>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_float32fe_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<float, true>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  // int 8 is not endian specific
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_int8_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<int8_t, true>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
  FAST_AUDIO_DATA_EXPORT int fad_render_sdl_from_uint8_frames(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb) {
    int res = convert_frames_to_doubles_swapped<uint8_t, true>(framebuffer);
    if (res) {
      return res;
    }
    TRY_BLOCK_START
      render_sdl_fft_data(pixels, pitch, height, sample_rate_hz, sr, sg, sb, er, eg, eb);
    TRY_BLOCK_END
    return res;
  }
}