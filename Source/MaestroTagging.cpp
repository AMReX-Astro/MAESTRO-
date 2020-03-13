
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// 
void
Maestro::RetagArray(const Box& bx,
                    const int lev, IntVector& tag_array)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RetagArray()", RetagArray);

    // re-compute tag_array since the actual grid structure changed due to buffering
    // this is required in order to compute numdisjointchunks, r_start_coord, r_end_coord

    Abort("Create a local copy of MaestroTagging.cpp in your build directory.\nHere is a sample that tags the refined region including buffer zone");

   // Tag on regions including buffer cells
   auto lo = bx.loVect3d()[AMREX_SPACEDIM-1];
   auto hi = bx.hiVect3d()[AMREX_SPACEDIM-1];
   const auto max_lev = max_radial_level + 1;

   for (auto r = lo; r <= hi; ++r) {
       tag_array[lev-1+max_lev*(r/2)] = TagBox::SET;
   }
}


void
Maestro::TagBoxes(TagBoxArray& tags, 
                  const MFIter& mfi,
                  const int lev, IntVector& tag_array,
                  const Real time)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TagBoxes()", TagBoxes);

    Abort("Create a local copy of MaestroTagging.cpp in your build directory.\nHere is a sample that tags the full grid using a predetermined tag array");

    const Array4<char> tag = tags.array(MFItInfo);
    const int * AMREX_RESTRICT tag_array_p = tag_array.dataPtr();
    const int max_lev = max_radial_level + 1;

    const Box& tilebox  = mfi.tilebox();
    const auto dx = geom[lev].CellSizeArray();

    AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        if (tag_array_p[lev+max_lev*r] > 0) {
            tag(i,j,k) = TagBox::SET;
        }
    });
}

// 
void
Maestro::StateError(TagBoxArray& tags, const MultiFab& state_mf, 
                   const MFIter& mfi,
                   const int lev, IntVector& tag_array,
                   const Real time)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::StateError()", StateError);

    Abort("Create a local copy of MaestroTagging.cpp in your build directory.\nHere is a sample that tags the temperature above 6.5e8");

    const Array4<char> tag = tags.array(mfi);
    const Array4<const Real> state = state_mf.array(mfi);
    int * AMREX_RESTRICT tag_array_p = tag_array.dataPtr();
    const int max_lev = max_radial_level + 1;

    const Box& tilebox  = mfi.tilebox();
    const auto dx = geom[lev].CellSizeArray();

    // Tag on regions of high temperature
    AMREX_PARALLEL_FOR_3D(tilebox, i, j, k, {
        if (state(i,j,k,Temp) >= 6.5e8) {
            int r = AMREX_SPACEDIM == 2 ? j : k;

            tag(i,j,k) = TagBox::SET;
            tag_array_p[lev+max_lev*r] = TagBox::SET;
        }
    });
}