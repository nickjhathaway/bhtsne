#pragma once
/*
 * MallocPtr.hpp
 *
 *  Created on: Mar 22, 2017
 *      Author: nick
 */



#include "bhtsne/common.h"


namespace bhtsne {

template <typename T>
class MallocPtr {
public:
  static std::shared_ptr<T> NumElements(const size_t numElements,
					const bool zero) {
    return NumBytes(sizeof(T) * numElements, zero);
  }

  static std::shared_ptr<T> NumBytes(const uint64_t numBytes,
				     const bool zero) {
    T* rawPtr = static_cast<T*>(malloc(numBytes));

    if (!rawPtr) {
      std::cerr << "could not allocate " << numBytes << " bytes of data";
      throw std::bad_alloc();
    }

    if(zero) {
      memset(rawPtr, 0, numBytes);
    }

    return std::shared_ptr<T>(rawPtr, [](T* ptr){ free(ptr); });
  }
};


}  // namespace bhtsne
