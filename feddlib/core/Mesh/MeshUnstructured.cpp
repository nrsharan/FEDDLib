#include "MeshUnstructured_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "MeshUnstructured_def.hpp"
namespace FEDD {
    template class MeshUnstructured<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

