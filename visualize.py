import pyaudio
import threading
import pygame
import struct
import ctypes
import sys
import os
try:
    from queue import Queue, Empty
except ImportError:
    from Queue import Queue, Empty
from os.path import dirname, join, isfile, abspath

pa = pyaudio.PyAudio()

output_device_name_try_list = []
input_device_name_try_list = []
base_dir = dirname(abspath(__file__))
input_preference_filename = join(base_dir, "idev_preference.txt")
output_preference_filename = join(base_dir, "odev_preference.txt")
if isfile(input_preference_filename):
    with open(input_preference_filename, "r") as fl:
        for line in fl:
            line = line.strip("\r\n")
            if len(line) == 0 or line.startswith('#'):
                continue
            input_device_name_try_list.append(line.lower())
else:
    with open(input_preference_filename, "w") as fl:
        fl.write("""# input device preference list
#   the first input device with a name that contains the full line (including
#   spaces, but is case insensitive) is matched for each line (except lines
#   starting with '#') starting with the first line. If no lines are matched,
#   then the first input device from the system is used.""")
if isfile(output_preference_filename):
    with open(output_preference_filename, "r") as fl:
        for line in fl:
            line = line.strip("\r\n")
            if len(line) == 0 or line.startswith('#'):
                continue
            output_device_name_try_list.append(line.lower())
else:
    with open(output_preference_filename, "w") as fl:
        fl.write("""# output device preference list
#   the first output device with a name that contains the full line (including
#   spaces, but is case insensitive) is matched for each line (except lines
#   starting with '#') starting with the first line. If no lines are matched,
#   then the first output device from the system is used.""")

def get_devices():
    devices = [None] * pa.get_device_count()
    for i in range(len(devices)):
        devices[i] = pa.get_device_info_by_index(i)
    return [
        (i, dev)
        for i, dev in enumerate(devices)
        if "@System32" not in dev["name"] and not dev["name"].endswith(" ()")
    ]

def open_stream_try_list(typ, try_list=None, verbose=False, **kwargs):
    assert typ in ["input", "output"]
    if verbose:
        print("trying to get an", typ, "device")
    if try_list is None:
        try_list = input_device_name_try_list if typ == "input" else output_device_name_try_list
    devices = get_devices()
    if typ == "output":
        kwargs["output"] = True
        narrow_devices = [(i, dev) for i, dev in devices if dev["maxInputChannels"] == 0 and dev["maxOutputChannels"] > 0]
    else:
        narrow_devices = [(i, dev) for i, dev in devices if dev["maxInputChannels"] > 0 and dev["maxOutputChannels"] == 0]
        kwargs["input"] = True
    if verbose:
        print("narrowed devices:\n  " + "\n  ".join(map(lambda x: "%u: %s" % (x[0], x[1]["name"]), narrow_devices)))
    for name in try_list:
        if verbose:
            print(f"using try_list entry '{name}' and searching over narrowed devices")
        found = [
            (i, dev)
            for i, dev in narrow_devices
            if name in dev["name"].lower()
        ]
        if verbose:
            print(f"  found {len(found)} results from searching")
        if len(found):
            if verbose:
                print("  found results:\n    " + "\n    ".join(map(lambda x: "%u: %s" % (x[0], x[1]["name"]), found)))
            print("using '%s' for '%s'" % (found[0][1]["name"], typ))
            kwargs[typ + "_device_index"] = found[0][0]
            break
    else:
        raise IOError("Could not find %s device that matches try_list %r" % (typ, try_list))
    return pa.open(**kwargs)



ostream = open_stream_try_list("output", rate=44100, channels=2, format=pyaudio.paInt16)
istream = open_stream_try_list("input", rate=44100, channels=2, format=pyaudio.paInt16)

def read_thread(numframes, q, istream, dct_sentinel):
    while dct_sentinel["running"]:
        q.put(istream.read(numframes))

def pygame_read_thread(numframes, istream, dct_sentinel):
    post = pygame.event.post
    Event = pygame.event.Event
    evt_type = pygame.USEREVENT
    while dct_sentinel["running"]:
        post(
            Event(
                evt_type,
                {
                    "framedata": istream.read(numframes)
                }
            )
        )

def write_thread(q, ostream, dct_sentinel, timeout=10):
    while dct_sentinel["running"]:
        try:
            ostream.write(q.get(timeout=timeout))
        except:
            pass

def passthru(istream, ostream, numframes=4410):
    dct_sentinel = {
        "running": True
    }
    q = Queue()
    t_read = threading.Thread(target=read_thread, args=(numframes, q, istream, dct_sentinel))
    t_write = threading.Thread(target=write_thread, args=(q, ostream, dct_sentinel))
    t_read.start()
    t_write.start()
    return [dct_sentinel, t_read, t_write, q]


def analyze(istream, ostream, numframes=4096,nchannels=2, sampwidth=2):
    dct_sentinel = {
        "running": True
    }
    iq = Queue()
    oq = Queue()
    typecode = {1: "b", 2: "h", 4: "i"}[sampwidth]
    s = struct.Struct(f"<{numframes * nchannels}{typecode}")
    t_read = threading.Thread(target=read_thread, args=(numframes, iq, istream, dct_sentinel))
    t_write = threading.Thread(target=write_thread, args=(oq, ostream, dct_sentinel))
    t_read.start()
    t_write.start()
    session = [dct_sentinel, t_read, t_write, iq, oq]
    pygame.display.init()
    width, height = 640, 480
    surf = pygame.display.set_mode((width, height))
    bkgr_color = (0, 0, 0)
    while dct_sentinel["running"]:
        events = pygame.event.get()
        for evt in events:
            if evt.type == pygame.QUIT:
                dct_sentinel["running"] = False
            elif evt.type == pygame.VIDEORESIZE:
                surf = pygame.display.set_mode(evt.size, pygame.RESIZABLE)
                width, height = evt.size
                surf.fill(bkgr_color)
        try:
            data = iq.get(block=False)
        except Empty:
            continue
        data1 = s.unpack(data)
        data2 = [0.0] * numframes
        for idx2, i in enumerate(range(0, len(data1), 2)):
            data2[idx2] = (data1[i] + data1[i + 1]) / 256 ** sampwidth
    t_read.join()
    t_write.join()


def init_fad():
    lib_path = join(base_dir, "win32/FastAudioData.dll" if os.name == "nt" else "linux/libFastAudioData.so")
    fad = ctypes.CDLL(lib_path)
    fad.fad_get_audio_buffer_size.restype = ctypes.c_size_t
    fad.fad_get_audio_buffer_size.argtypes = []
    endian = "fe" if sys.byteorder == "big" else "ne"
    fn_mappings = {
        pyaudio.paInt32: f"fad_render_sdl_from_int32{endian}_frames",
        pyaudio.paInt16: f"fad_render_sdl_from_int16{endian}_frames",
        pyaudio.paFloat32: f"fad_render_sdl_from_float32{endian}_frames",
        pyaudio.paInt8: "fad_render_sdl_from_int8_frames",
        pyaudio.paUInt8: "fad_render_sdl_from_uint8_frames"
    }
    fn_names = [
        f"fad_render_sdl_from_int32{endian}_frames",
        f"fad_render_sdl_from_int16{endian}_frames",
        f"fad_render_sdl_from_uint32{endian}_frames",
        f"fad_render_sdl_from_uint16{endian}_frames",
        f"fad_render_sdl_from_float64{endian}_frames",
        f"fad_render_sdl_from_float32{endian}_frames",
        "fad_render_sdl_from_int8_frames",
        "fad_render_sdl_from_uint8_frames"
    ]
    # c_uint8_p = ctypes.POINTER(ctypes.c_uint8)
    # ptr_type = ctypes.c_void_p
    ptr_type = ctypes.c_char_p
    for k in fn_names:
        fn = getattr(fad, k)
        fn.restype = ctypes.c_int
        fn.argtypes = [ptr_type, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_uint8, ctypes.c_uint8, ctypes.c_uint8, ctypes.c_uint8, ctypes.c_uint8, ctypes.c_uint8]
        fn.__doc__ = f"void {k}(const uint8_t *framebuffer, void *pixels, size_t pitch, size_t height, size_t sample_rate_hz, uint8_t sr, uint8_t sg, uint8_t sb, uint8_t er, uint8_t eg, uint8_t eb)"
    fad.fad_build_buffers.restype = None
    fad.fad_build_buffers.argtypes = [ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t]
    fad.fad_build_buffers.__doc__ = "void fad_build_buffers(size_t framesize, size_t nchannels, size_t fft_size, size_t display_width)"
    fad.fad_destroy_buffers.restype = None
    fad.fad_destroy_buffers.argtypes = []
    fad.fad_destroy_buffers.__doc__ = "void fad_destroy_buffers()"
    fad.fad_set_averaging.restype = None
    fad.fad_set_averaging.argtypes = [ctypes.c_size_t]
    fad.fad_set_averaging.__doc__ = "void fad_set_averaging(size_t averaging)"
    fad.fad_get_averaging.restype = ctypes.c_size_t
    fad.fad_get_averaging.argtypes = []
    fad.fad_get_averaging.__doc__ = "size_t fad_get_averaging()"
    return fad, fn_mappings


def analyze_pygame(istream, ostream, numframes=2048, nchannels=2, sampwidth=2):
    dct_sentinel = {
        "running": True
    }
    assert istream._rate == ostream._rate
    assert istream._format == ostream._format
    sample_rate_hz = istream._rate
    oq = Queue()
    t_read = threading.Thread(target=pygame_read_thread, args=(numframes, istream, dct_sentinel))
    t_write = threading.Thread(target=write_thread, args=(oq, ostream, dct_sentinel, 1))
    pygame.display.init()
    width, height = 640, 480
    window_flags = pygame.HWSURFACE|pygame.DOUBLEBUF|pygame.RESIZABLE
    # window_flags = pygame.RESIZABLE
    surf = pygame.display.set_mode((width, height), window_flags)
    assert surf.get_masks()[:3] == (0xFF0000, 0xFF00, 0xFF)
    bkgr_color = (0, 0, 0)
    t_read.start()
    t_write.start()
    fad, fn_mappings = init_fad()
    fad.fad_build_buffers(sampwidth, nchannels, numframes >> 1, width)
    fad_render_sdl_from_frames = getattr(fad, fn_mappings[istream._format])
    fad_audio_bufsize = fad.fad_get_audio_buffer_size()
    c_uint8_p = ctypes.POINTER(ctypes.c_uint8)
    is_paused = False
    while dct_sentinel["running"]:
        evt = pygame.event.wait()
        events = pygame.event.get()
        events.insert(0, evt)
        del evt
        audio_events_encountered = 0
        for evt in events:
            if evt.type == pygame.QUIT:
                dct_sentinel["running"] = False
            elif evt.type == pygame.VIDEORESIZE:
                surf = pygame.display.set_mode(evt.size, window_flags)
                assert surf.get_masks()[:3] == (0xFF0000, 0xFF00, 0xFF)
                width, height = evt.size
                surf.fill(bkgr_color)
                fad.fad_build_buffers(sampwidth, nchannels, numframes >> 1, width)
                fad_audio_bufsize = fad.fad_get_audio_buffer_size()
            elif evt.type == pygame.USEREVENT:
                data = evt.framedata
                if audio_events_encountered > 0:
                    print("more than one audio event encountered in the same event cycle, dropping")
                    continue
                oq.put(data)
                audio_events_encountered += 1
                assert len(data) == fad_audio_bufsize
                # sys.stdout.write('pyth hex (first 108 bytes):')
                # for i in range(108):
                #     sys.stdout.write(f" {digits[data[i] >> 4]}{digits[data[i] & 0xF]}")
                # sys.stdout.write('\n')
                if not is_paused:
                    surf.lock()
                    try:
                        res = fad_render_sdl_from_frames(data, surf._pixels_address, surf.get_pitch(), height, sample_rate_hz, 0, 0, 255, 0, 255, 0)
                    except:
                        print("pixels_address = 0x%016X" % surf._pixels_address)
                        raise
                    finally:
                        surf.unlock()
                    if res:
                        print("WARN: error code", res)
                    pygame.display.flip()
            elif evt.type == pygame.KEYDOWN:
                if evt.key == pygame.K_UP:
                    cur_avg = fad.fad_get_averaging()
                    if cur_avg < 10:
                        print(f"Increasing averaging from {cur_avg} to {cur_avg + 1}")
                        fad.fad_set_averaging(cur_avg + 1)
                elif evt.key == pygame.K_DOWN:
                    cur_avg = fad.fad_get_averaging()
                    if cur_avg > 0:
                        print(f"Decreasing averaging from {cur_avg} to {cur_avg - 1}")
                        fad.fad_set_averaging(cur_avg - 1)
                elif evt.key == pygame.K_r:
                    print(f"Resetting averaging to 0")
                    fad.fad_set_averaging(0)
                elif evt.key == pygame.K_SPACE:
                    is_paused = not is_paused
    print("waiting for read thread to join")
    t_read.join()
    print("waiting for write thread to join")
    t_write.join()
    print("deallocating buffers")
    fad.fad_destroy_buffers()
    print("quiting pygame")
    pygame.quit()


if __name__ == "__main__":
    analyze_pygame(istream, ostream)