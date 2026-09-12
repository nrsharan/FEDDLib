#include "AssembleFE_NeoHooke2D_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFE_NeoHooke2D_def.hpp"
namespace FEDD {
    template class AssembleFE_NeoHooke2D<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
