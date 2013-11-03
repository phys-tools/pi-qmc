#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "Permutation.h"

Permutation::Permutation(const int n) : permutation(n) {
  for (int i=0; i<n; ++i) permutation(i)=i;
}

Permutation::Permutation(const Permutation& p)
  : permutation(p.permutation.copy()) {
}

Permutation& Permutation::operator=(const Permutation& p) {
  permutation=p.permutation;
  return *this;
}

void Permutation::reset() {
  for (int i=0; i<permutation.size(); ++i) permutation(i)=i;
}

Permutation& Permutation::prepend(const Permutation& p) {
  int *temp;
  temp = new int[permutation.size()];
  for (int i=0; i<permutation.size(); ++i) {
    temp[i]=permutation(p.permutation(i));
  }
  for (int i=0; i<permutation.size(); ++i) permutation(i)=temp[i];
  delete temp;
  return *this;
}

Permutation& Permutation::append(const Permutation& p) {
  int *temp;
  temp = new int[permutation.size()];
  for (int i=0; i<permutation.size(); ++i) {
    temp[i]=p.permutation(permutation(i));
  }
  for (int i=0; i<permutation.size(); ++i) permutation(i)=temp[i];
  delete temp;
  return *this;
}

void Permutation::setToInverse(const Permutation& p) {
  for (int i=0; i<permutation.size(); ++i) {
    permutation(p.permutation(i)) = i;
  }
}

bool Permutation::isIdentity() const {
  for (int i=0; i<permutation.size(); ++i) {
    if (permutation(i)!=i) return false;
  }
  return true;
}

int Permutation::inversionCount() const {
    // Note: should be order(n log n), but this is order(n^2)
    int count = 0;
    for (int i = 0; i < permutation.size() - 1; ++i) {
        for (int j = i + 1; j < permutation.size(); ++j) {
            if (permutation(i) > permutation(j)) {
                ++count;
            }
        }
    }
    return count;
}

int Permutation::sign() const {
    return ((inversionCount() & 1) == 0) ? 1 : -1;
}

std::ostream& operator<<(std::ostream &os, const Permutation &p) {
  return os << p.permutation;
}
