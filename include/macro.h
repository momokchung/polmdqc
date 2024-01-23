// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "darray.h"
#include "precision.h"
#include <omp.h>

#ifndef MDQC_EXTERN_DEFINITION_FILE
#define MDQC_EXTERN_DEFINITION_FILE 0
#endif
#if MDQC_EXTERN_DEFINITION_FILE
#define MDQC_EXTERN
#else
#define MDQC_EXTERN extern
#endif

#define MAYBE_UNUSED [[maybe_unused]]

#if defined(__INTEL_COMPILER)
#define TINKER_ICPC
#endif
