
#ifndef SOME_LIB_EXPORT_H
#define SOME_LIB_EXPORT_H

#ifdef SOME_LIB_STATIC_DEFINE
#  define SOME_LIB_EXPORT
#  define SOME_LIB_NO_EXPORT
#else
#  ifndef SOME_LIB_EXPORT
#    ifdef some_lib_EXPORTS
        /* We are building this library */
#      define SOME_LIB_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define SOME_LIB_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef SOME_LIB_NO_EXPORT
#    define SOME_LIB_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef SOME_LIB_DEPRECATED
#  define SOME_LIB_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef SOME_LIB_DEPRECATED_EXPORT
#  define SOME_LIB_DEPRECATED_EXPORT SOME_LIB_EXPORT SOME_LIB_DEPRECATED
#endif

#ifndef SOME_LIB_DEPRECATED_NO_EXPORT
#  define SOME_LIB_DEPRECATED_NO_EXPORT SOME_LIB_NO_EXPORT SOME_LIB_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef SOME_LIB_NO_DEPRECATED
#    define SOME_LIB_NO_DEPRECATED
#  endif
#endif

#endif /* SOME_LIB_EXPORT_H */
