#ifndef BaseState_H_
#define BaseState_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>

template <class T>
class BaseState;

template <typename T>
struct BaseStateArray {
    T* AMREX_RESTRICT dptr;

    int nlev;
    int len;
    int nvar;

    /// default constructor
    AMREX_GPU_HOST_DEVICE
    constexpr BaseStateArray() noexcept
        : dptr(nullptr), nlev(0), len(0), nvar(0) {}

    /// copy constructor
    AMREX_GPU_HOST_DEVICE
    constexpr BaseStateArray(BaseStateArray<T> const& rhs) noexcept
        : dptr(rhs.dptr), nlev(rhs.nlev), len(rhs.len), nvar(rhs.nvar) {}

    /// initialize from pointer
    AMREX_GPU_HOST_DEVICE
    constexpr BaseStateArray(T* ptr, const int num_levs = 1,
                             const int length = 1, const int ncomp = 1) noexcept
        : dptr(ptr), nlev(num_levs), len(length), nvar(ncomp) {}

    AMREX_GPU_HOST_DEVICE
    void init(T* ptr, const int num_levs = 1, const int length = 1,
              const int ncomp = 1) noexcept {
        dptr = ptr;
        nlev = num_levs;
        len = length;
        nvar = ncomp;
    }

    void init(BaseState<T>& base_state) noexcept {
        dptr = base_state.dataPtr();
        nlev = base_state.nLevels();
        len = base_state.length();
        nvar = base_state.nComp();
    }

    AMREX_GPU_HOST_DEVICE
    explicit operator bool() const noexcept { return dptr != nullptr; }

    /// access elements. If only 2 arguments are passed in, component is set to 0
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& operator()(
        const int lev, const int i, const int n = 0) const noexcept {
        AMREX_ASSERT(lev >= 0);
        AMREX_ASSERT(i < this->len && i >= 0);
        AMREX_ASSERT(n < this->nvar && n >= 0);

        return dptr[(lev * len + i) * nvar + n];
    }

    /// if access using only one index, then this ignores the underlying data structure. Useful for e.g. applying the same operation to all data.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& operator()(
        const int i) const noexcept {
        return dptr[i];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& operator[](int i) noexcept {
        AMREX_ASSERT(i >= 0 && i < nlev * len * nvar);

        return dptr[i];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& operator[](
        const int i) const noexcept {
        AMREX_ASSERT(i >= 0 && i < nlev * len * nvar);

        return dptr[i];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T* ptr(int lev, int i = 0,
                                                    int n = 0) const noexcept {
        return dptr + (lev * len + i) * nvar + n;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T* dataPtr() noexcept {
        return this->dptr;
    }

    /// number of levels in base
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int nLevels() const noexcept {
        return nlev;
    }

    /// number of cells in base
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int length() const noexcept {
        return len;
    }

    /// number of componenets
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int nComp() const noexcept {
        return nvar;
    }
};

/// create a BaseStateArray
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE BaseStateArray<T> makeBaseStateArray(
    T* dptr, const int num_levs = 1, const int length = 1,
    const int ncomp = 1) noexcept {
    return BaseStateArray<T>{dptr, num_levs, length, ncomp};
}

template <class T>
class BaseState {
   public:
    /*
      public member functions
    */

    BaseState() : nlev(0), len(0), nvar(0) {}

    explicit BaseState(const int num_levs, const int length = 1,
                       const int ncomp = 1, const T val = T());

    /// create a BaseState object from an existing Fortran-friendly base state
    BaseState(const amrex::Gpu::ManagedVector<T>& src, const int num_levs,
              const int length = 1, const int ncomp = 1);

    BaseState(const amrex::Vector<T>& src, const int num_levs,
              const int length = 1, const int ncomp = 1);

    /// copy constructor. This makes a deep copy of the src.
    BaseState(const BaseState<T>& src);

    /// return a BaseStateArray object to allow for accessing
    /// the underlying data
    AMREX_FORCE_INLINE
    BaseStateArray<T const> array() const noexcept {
        return makeBaseStateArray<T const>(base_data.dataPtr(), nlev, len,
                                           nvar);
    }

    /// return a BaseStateArray object to allow for accessing
    /// the underlying data
    AMREX_FORCE_INLINE
    BaseStateArray<T> array() noexcept {
        return makeBaseStateArray<T>(base_data.dataPtr(), nlev, len, nvar);
    }

    AMREX_FORCE_INLINE
    BaseStateArray<const T> const_array() const noexcept {
        return makeBaseStateArray<const T>(base_data.dataPtr(), nlev, len,
                                           nvar);
    }

    /// return ith element of underlying data array
    AMREX_FORCE_INLINE
    T array(const int i) noexcept {
        BaseStateArray<T> arr =
            makeBaseStateArray<T>(base_data.dataPtr(), nlev, len, nvar);
        return arr(i);
    }

    /// allocate memory for the BaseState
    void define(const int num_levs, const int length = 1, const int ncomp = 1,
                const T val = T());

    void resize(const int num_levs, const int length = 1, const int ncomp = 1,
                const T val = T());

    /// set to some scalar value
    void setVal(const T& val);

    /// set the comp'th component to some scalar value
    void setVal(const int comp, const T& val);

    /// number of levels in base
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int nLevels() const noexcept {
        return nlev;
    };

    /// number of cells in base
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int length() const noexcept {
        return len;
    };

    /// number of components
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int nComp() const noexcept {
        return nvar;
    };

    /// deep copy
    void copy(const BaseState<T>& src);
    void copy(const BaseStateArray<T> src);
    void copy(const amrex::Gpu::ManagedVector<T>& src);
    void copy(const amrex::Vector<T>& src);

    void toVector(amrex::Vector<T>& vec) const;
    void toVector(amrex::Gpu::ManagedVector<T>& vec) const;

    /// swap the data with src
    void swap(BaseState<T>& src);

    T* dataPtr() noexcept { return base_data.dataPtr(); };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T* dataPtr() const noexcept {
        return this->dptr;
    }

    /// scalar addition to the whole base state
    template <class U>
    friend BaseState<U> operator+(const U val, const BaseState<U>& p);
    BaseState<T>& operator+=(const T val);

    /// element-wise addition
    template <class U>
    friend BaseState<U> operator+(const BaseState<U>& lhs,
                                  const BaseState<U>& rhs);
    BaseState<T>& operator+=(const BaseState<T>& rhs);

    /// scalar subtraction from the whole base state
    template <class U>
    friend BaseState<U> operator-(const U val, const BaseState<U>& p);
    BaseState<T>& operator-=(const T val);

    /// element-wise subtraction
    template <class U>
    friend BaseState<U> operator-(const BaseState<U>& lhs,
                                  const BaseState<U>& rhs);
    BaseState<T>& operator-=(const BaseState<T>& rhs);

    /// scalar multiplication of the whole base state
    template <class U>
    friend BaseState<U> operator*(const U val, const BaseState<U>& p);
    BaseState<T>& operator*=(const T val);

    /// element-wise multiplication
    template <class U>
    friend BaseState<U> operator*(const BaseState<U>& lhs,
                                  const BaseState<U>& rhs);
    BaseState<T>& operator*=(const BaseState<T>& rhs);

    /// scalar division of the whole base state
    template <class U>
    friend BaseState<U> operator/(const U val, const BaseState<U>& p);
    BaseState<T> operator/(const T val);
    BaseState<T>& operator/=(const T val);
    BaseState<T> operator/(const T val) const noexcept;

    /// element-wise division
    template <class U>
    friend BaseState<U> operator/(const BaseState<U>& lhs,
                                  const BaseState<U>& rhs);
    BaseState<T>& operator/=(const BaseState<T>& rhs);

    /// comparison operator
    template <class U>
    friend bool operator==(const BaseState<U>& lhs, const BaseState<U>& rhs);

    /// comparison operator
    template <class U>
    friend bool operator!=(const BaseState<U>& lhs, const BaseState<U>& rhs);

   protected:
    amrex::Gpu::ManagedVector<T> base_data;

    int nlev;
    int len;
    int nvar;
};

template <class T>
BaseState<T>::BaseState(const int num_levs, const int length, const int ncomp,
                        const T val)
    : nlev(num_levs), len(length), nvar(ncomp) {
    base_data.resize(nlev * len * nvar);
    base_data.shrink_to_fit();
    std::fill(base_data.begin(), base_data.end(), val);
}

/// make sure that this class is never instantiated with bool
template <>
BaseState<bool>::BaseState(const int, const int, const int,
                           const bool) = delete;

template <class T>
BaseState<T>::BaseState(const amrex::Gpu::ManagedVector<T>& src,
                        const int num_levs, const int length, const int ncomp)
    : nlev(num_levs), len(length), nvar(ncomp) {
    base_data.resize(nlev * len * nvar);
    base_data.shrink_to_fit();
    // need to switch from Fortran array-ordering to C++ ordering
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto n = 0; n < nvar; ++n) {
                base_data[n + nvar * (r + len * l)] =
                    src[l + nlev * (r + len * n)];
            }
        }
    }
}

template <class T>
BaseState<T>::BaseState(const amrex::Vector<T>& src, const int num_levs,
                        const int length, const int ncomp)
    : nlev(num_levs), len(length), nvar(ncomp) {
    base_data.resize(nlev * len * nvar);
    base_data.shrink_to_fit();
    // need to switch from Fortran array-ordering to C++ ordering
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto n = 0; n < nvar; ++n) {
                base_data[n + nvar * (r + len * l)] =
                    src[l + nlev * (r + len * n)];
            }
        }
    }
}

template <class T>
BaseState<T>::BaseState(const BaseState<T>& src)
    : nlev(src.nLevels()), len(src.length()), nvar(src.nComp()) {
    base_data.resize(nlev * len * nvar);
    base_data.shrink_to_fit();
    BaseStateArray<T> base_arr = this->array();
    const BaseStateArray<const T> src_arr = src.array();
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto n = 0; n < nvar; ++n) {
                base_arr(l, r, n) = src_arr(l, r, n);
            }
        }
    }
}

template <class T>
void BaseState<T>::define(const int num_levs, const int length, const int ncomp,
                          const T val) {
    nlev = num_levs;
    len = length;
    nvar = ncomp;

    base_data.resize(nlev * length * nvar);
    base_data.shrink_to_fit();
    std::fill(base_data.begin(), base_data.end(), val);
}

template <class T>
void BaseState<T>::resize(const int num_levs, const int length, const int ncomp,
                          const T val) {
    if (base_data.size() > 0) {
        base_data.resize(0);
    }
    this->define(num_levs, length, ncomp);
}

template <class T>
void BaseState<T>::setVal(const T& val) {
    BaseStateArray<T> base_arr = this->array();
    AMREX_PARALLEL_FOR_1D(nvar * len * nlev, i, { base_arr(i) = val; });

    amrex::Gpu::synchronize();
}

template <class T>
void BaseState<T>::setVal(const int comp, const T& val) {
    BaseStateArray<T> base_arr = this->array();
    AMREX_PARALLEL_FOR_1D(len * nlev, i, { base_arr(0, i, comp) = val; });

    amrex::Gpu::synchronize();
}

template <class T>
void BaseState<T>::copy(const BaseState<T>& src) {
    AMREX_ASSERT(nlev == src.nlev);
    AMREX_ASSERT(nvar == src.nvar);
    AMREX_ASSERT(len == src.len);

    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto comp = 0; comp < nvar; ++comp) {
                base_data[comp + nvar * (r + len * l)] =
                    src.base_data[comp + nvar * (r + len * l)];
            }
        }
    }
}

template <class T>
void BaseState<T>::copy(const BaseStateArray<T> src) {
    AMREX_ASSERT(nlev == src.nlev);
    AMREX_ASSERT(nvar == src.nvar);
    AMREX_ASSERT(len == src.len);

    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto comp = 0; comp < nvar; ++comp) {
                base_data[comp + nvar * (r + len * l)] = src(l, r, comp);
            }
        }
    }
}

template <class T>
void BaseState<T>::copy(const amrex::Gpu::ManagedVector<T>& src) {
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto comp = 0; comp < nvar; ++comp) {
                base_data[comp + nvar * (r + len * l)] =
                    src[l + nlev * (r + len * comp)];
            }
        }
    }
}

template <class T>
void BaseState<T>::copy(const amrex::Vector<T>& src) {
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto comp = 0; comp < nvar; ++comp) {
                base_data[comp + nvar * (r + len * l)] =
                    src[l + nlev * (r + len * comp)];
            }
        }
    }
}

template <class T>
void BaseState<T>::toVector(amrex::Vector<T>& vec) const {
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto comp = 0; comp < nvar; ++comp) {
                vec[l + nlev * (r + len * comp)] =
                    base_data[comp + nvar * (r + len * l)];
            }
        }
    }
}

template <class T>
void BaseState<T>::toVector(amrex::Gpu::ManagedVector<T>& vec) const {
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto comp = 0; comp < nvar; ++comp) {
                vec[l + nlev * (r + len * comp)] =
                    base_data[comp + nvar * (r + len * l)];
            }
        }
    }
}

template <class T>
void BaseState<T>::swap(BaseState<T>& src) {
    AMREX_ASSERT(nlev == src.nlev);
    AMREX_ASSERT(nvar == src.nvar);
    AMREX_ASSERT(len == src.len);

    auto arr = this->array();
    auto src_arr = src.array();

    AMREX_PARALLEL_FOR_1D(len * nlev * nvar, i, {
        T temp = arr(i);
        arr(i) = src_arr(i);
        src_arr(i) = temp;
    });

    amrex::Gpu::synchronize();
}

template <class T>
BaseState<T> operator+(const T val, const BaseState<T>& p) {
    // make a deep copy
    BaseState<T> s(p);
    s += val;
    return s;
}

template <class T>
BaseState<T>& BaseState<T>::operator+=(const T val) {
    BaseStateArray<T> base_arr = this->array();
    AMREX_PARALLEL_FOR_1D(nvar * len * nlev, i, { base_arr(i) += val; });
    amrex::Gpu::synchronize();
    return *this;
}

template <class T>
BaseState<T> operator+(const BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    BaseState<T> s(lhs);
    s += rhs;
    return s;
}

template <class T>
BaseState<T>& BaseState<T>::operator+=(const BaseState<T>& rhs) {
    AMREX_ASSERT(nlev == rhs.nlev);
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    BaseStateArray<T> base_arr = this->array();
    const BaseStateArray<const T> rhs_arr = rhs.array();
    AMREX_PARALLEL_FOR_1D(nvar * len * nlev, i, { base_arr(i) += rhs_arr(i); });
    amrex::Gpu::synchronize();
    return *this;
}

template <class T>
BaseState<T> operator-(const T val, const BaseState<T>& p) {
    // make a deep copy
    BaseState<T> s(p);
    s -= val;
    return s;
}

template <class T>
BaseState<T>& BaseState<T>::operator-=(const T val) {
    BaseStateArray<T> base_arr = this->array();
    AMREX_PARALLEL_FOR_1D(nvar * len * nlev, i, { base_arr(i) -= val; });
    amrex::Gpu::synchronize();
    return *this;
}

template <class T>
BaseState<T> operator-(const BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    BaseState<T> s(lhs);
    s -= rhs;
    return s;
}

template <class T>
BaseState<T>& BaseState<T>::operator-=(const BaseState<T>& rhs) {
    AMREX_ASSERT(nlev == rhs.nlev);
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    BaseStateArray<T> base_arr = this->array();
    const BaseStateArray<const T> rhs_arr = rhs.array();
    AMREX_PARALLEL_FOR_1D(nvar * len * nlev, i, { base_arr(i) -= rhs_arr(i); });
    amrex::Gpu::synchronize();
    return *this;
}

template <class T>
BaseState<T> operator*(const T val, const BaseState<T>& p) {
    // make a deep copy
    BaseState<T> s(p);
    s *= val;
    return s;
}

template <class T>
BaseState<T>& BaseState<T>::operator*=(const T val) {
    BaseStateArray<T> base_arr = this->array();
    AMREX_PARALLEL_FOR_1D(nvar * len * nlev, i, { base_arr(i) *= val; });
    amrex::Gpu::synchronize();
    return *this;
}

template <class T>
BaseState<T> operator*(const BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    BaseState<T> s(lhs);
    s *= rhs;
    return s;
}

template <class T>
BaseState<T>& BaseState<T>::operator*=(const BaseState<T>& rhs) {
    AMREX_ASSERT(nlev == rhs.nlev);
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    BaseStateArray<T> base_arr = this->array();
    const BaseStateArray<const T> rhs_arr = rhs.array();
    AMREX_PARALLEL_FOR_1D(nvar * len * nlev, i, { base_arr(i) *= rhs_arr(i); });
    amrex::Gpu::synchronize();
    return *this;
}

template <class T>
BaseState<T> operator/(const T val, const BaseState<T>& p) {
    // make a deep copy
    BaseState<T> s(p);
    s /= val;
    return s;
}

template <class T>
BaseState<T> BaseState<T>::operator/(const T val) {
    // make a deep copy
    BaseState<T> s(*this);
    s /= val;
    return s;
}

template <class T>
BaseState<T>& BaseState<T>::operator/=(const T val) {
    BaseStateArray<T> base_arr = this->array();
    AMREX_PARALLEL_FOR_1D(nvar * len * nlev, i, { base_arr(i) /= val; });
    amrex::Gpu::synchronize();
    return *this;
}

template <class T>
BaseState<T> BaseState<T>::operator/(const T val) const noexcept {
    BaseState<T> s(nlev, len, nvar);
    BaseStateArray<T> s_arr = s.array();
    BaseStateArray<T> base_arr = this->array();
    AMREX_PARALLEL_FOR_1D(nvar * len * nlev, i,
                          { s_arr(i) = base_arr(i) / val; });
    amrex::Gpu::synchronize();
    return s;
}

template <class T>
BaseState<T> operator/(const BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    BaseState<T> s(lhs);
    s /= rhs;
    return s;
}

template <class T>
BaseState<T>& BaseState<T>::operator/=(const BaseState<T>& rhs) {
    AMREX_ASSERT(nlev == rhs.nlev);
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    BaseStateArray<T> base_arr = this->array();
    const BaseStateArray<const T> rhs_arr = rhs.array();
    AMREX_PARALLEL_FOR_1D(nvar * len * nlev, i, { base_arr(i) /= rhs_arr(i); });
    amrex::Gpu::synchronize();
    return *this;
}

template <class T>
bool operator==(const BaseState<T>& lhs, const BaseState<T>& rhs) {
    AMREX_ASSERT(lhs.nlev == rhs.nlev);
    AMREX_ASSERT(lhs.nvar == rhs.nvar);
    AMREX_ASSERT(lhs.len == rhs.len);

    const BaseStateArray<const T> lhs_arr = lhs.array();
    const BaseStateArray<const T> rhs_arr = rhs.array();

    for (auto i = 0; i < lhs.nvar * lhs.len * lhs.nlev; ++i) {
        if (lhs_arr(i) != rhs_arr(i)) {
            return false;
        }
    }
    return true;
}

template <class T>
bool operator!=(const BaseState<T>& lhs, const BaseState<T>& rhs) {
    AMREX_ASSERT(lhs.nlev == rhs.nlev);
    AMREX_ASSERT(lhs.nvar == rhs.nvar);
    AMREX_ASSERT(lhs.len == rhs.len);

    const BaseStateArray<const T> lhs_arr = lhs.array();
    const BaseStateArray<const T> rhs_arr = rhs.array();

    for (auto i = 0; i < lhs.nvar * lhs.len * lhs.nlev; ++i) {
        if (lhs_arr(i) != rhs_arr(i)) {
            return true;
        }
    }
    return false;
}
#endif
