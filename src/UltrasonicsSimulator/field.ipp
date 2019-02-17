


template<class T>
Field<T>::Field(const std::array<size_t, 3> N, const std::array<double, 3> L, const std::array<double, 3> origin):
    _N(N), _L(L), orig(origin)
{

    values.resize(_N[0] * _N[1] * _N[2]);
    
}


template<class T>
const T& Field<T>::pos(const std::array<size_t, 3> c) const {

    const size_t x = c[0];
    const size_t y = c[1];
    const size_t z = c[2];

    const size_t index = x + _N[0] * y + _N[0] * _N[1] * z;

    return values[index];

}



template<class T>
T& Field<T>::pos(const std::array<size_t, 3> c) {

    return const_cast<T &>(static_cast<const Field<T> &>(*this).pos(c));

}



template<class T>
const T& Field<T>::operator() (const std::array<size_t, 3> c) const {

    return this->pos(c);

}



template<class T>
T& Field<T>::operator() (const std::array<size_t, 3> c) {

    return this->pos(c);

}



template<class T>
const std::array<size_t, 3> Field<T>::N() const {
    return _N;
}



template<class T>
const size_t Field<T>::N(const size_t i) const {
    return _N[i];
}



template<class T>
const std::array<double, 3> Field<T>::L() const {
    return _L;
}



template<class T>
const double Field<T>::L(const size_t i) const {
    return _L[i];
}
