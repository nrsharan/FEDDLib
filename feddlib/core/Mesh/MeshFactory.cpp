#include "MeshFactory_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "MeshFactory_def.hpp"
namespace FEDD {
    template class MeshFactory<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
