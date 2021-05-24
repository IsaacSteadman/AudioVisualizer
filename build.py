import os
from subprocess import call
from sys import argv
import json
from os.path import join, dirname, splitext, abspath
from typing import List


def p_call(args, **kwargs):
    print(args, dict(kwargs))
    call(args, **kwargs)


base_dir = abspath(dirname(__file__))
obj_name = "FastAudioData"


def get_sources(d: str) -> List[str]:
    print("get_sources(%r)" % d)
    res = []
    for f in os.listdir(d):
        base, ext = splitext(f)
        if ext in [".c", ".cpp"]:
            print("  " + f)
            res.append(join(d, f))
    return res


lst_sources = get_sources(base_dir)


class Config(object):
    def __init__(self, json: dict):
        self.json = json
    
    def resolve(self, k: str) -> str:
        s = self.json[k]
        return self.resolve_internal(s)
    
    def resolve_internal(self, s: str) -> str:
        pos1 = s.find("<")
        pos2 = s.find(">", pos1)
        while pos1 >= 0 and pos2 >= 0:
            s = s[:pos1] + self.resolve(s[pos1 + 1:pos2]) + s[pos2 + 1:]
            pos1 = s.find("<")
            pos2 = s.find(">", pos1 + 1)
        return s
    
    def resolve_dict_list_str_reduce(self, k: str, reduce_fn):
        dct: dict = self.json[k]
        assert isinstance(dct, dict)
        return {k1: reduce_fn([self.resolve_internal(v1) for v1 in v]) for k1, v in dct.items()}


if os.name == "nt":
    with open(join(base_dir, "win32", "build_config.json"), "r") as fl:
        config = Config(json.load(fl))
    msvc_dir = config.resolve("msvc_dir")
    win_kit_dir = config.resolve("win_kit_dir")
    compiler_args = [
        config.resolve("cl.exe"), "/Wall",
        "/D", "NDEBUG", "/D", "_WINDOWS", "/D", "_USRDLL",
        "/D", "_WINDLL", "/D", "_UNICODE", "/D", "UNICODE",
        "/std:c++17", "/DEBUG"
    ]
    dylib_args = [
        "/LD",
    ]
    at_end = [
        "/link", "/MACHINE:X64", "/DEBUG",
        "/OUT:%s" % join(base_dir, "win32", "./%s.dll" % obj_name),
        "/SUBSYSTEM:WINDOWS",
        "/DYNAMICBASE", "kernel32.lib", "user32.lib",
        "gdi32.lib", "winspool.lib", "comdlg32.lib",
        "advapi32.lib", "shell32.lib", "ole32.lib",
        "oleaut32.lib", "uuid.lib", "odbc32.lib", "odbccp32.lib"
    ]
    env = dict(os.environ, **config.resolve_dict_list_str_reduce("cl_env", lambda x: ";".join(x)))
elif os.name == "posix":
    compiler_args = [
        "gcc", "-fdiagnostics-color=always", "-lstdc++",
        "-std=c++14", "-ldl", "-pthread", "-L.", "-Wall"
    ]
    dylib_args = [
        "-shared", "-fPIC", "-o", join(base_dir, "linux", "lib%s.so" % obj_name)
    ]
    at_end = []
else:
    raise OSError("OS not recognized")

if len(lst_sources) == 0:
    print("No sources found")
    raise SystemExit(-1)

compiler_args.extend(dylib_args)
compiler_args.extend(lst_sources)
compiler_args.extend(at_end)
if __name__ == "__main__":
    if len(argv) > 1:
        if len(argv) == 2 and argv[1].lower() == "clean":
            if os.name == "nt":
                path = join(base_dir, "win32")
                for f in os.listdir(path):
                    base, ext = splitext(f)
                    if ext in [".exp", ".lib", ".obj", ".dll"]:
                        os.remove(join(path, f))
                        print("removing " + f)
            else:
                print("Unsupported operating system: '%s'" % os.name)
        else:
            print("bad args, only supports no args or 'clean'")
    else:
        # print("Compilargs\n  " + "\n  ".join(compiler_args))
        call(compiler_args, cwd=join(base_dir, "win32"), env=env)
