


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
const T Field<T>::interp(const std::array<double, 3> c) const {

    size_t x, y, z;
    
    if (this->N(0) <= 1) {
	x = orig[0] + this->L(0)/2;
    }
    else {
	x = std::floor(std::min(std::max((c[0] - orig[0]) / (this->L(0) / this->N(0)), 0.0), static_cast<double>(this->N(0)-1)));
    }
    if (this->N(1) <= 1) {
	y = orig[1] + this->L(1)/2;
    }
    else {
	y = std::floor(std::min(std::max((c[1] - orig[1]) / (this->L(1) / this->N(1)), 0.0), static_cast<double>(this->N(1)-1)));
    }
    if (this->N(2) <= 1) {
	z = orig[2] + this->L(2)/2;
    }
    else {
	z = std::floor(std::min(std::max((c[2] - orig[2]) / (this->L(2) / this->N(2)), 0.0), static_cast<double>(this->N(2)-1)));
    }

    //std::cout << "a " << (orig[1] + this->L(1) + c[1]) << " " <<  (this->L(1) / this->N(1)) << "\n";
    //std::cout << "x: " << x << ", y: " << y << ", z: " << z << " " << std::floor(static_cast<double>(this->N(2)-1)) << "\n";
    
    return this->pos({x, y, z});
    
}



template<class T>
std::array<size_t, 3> Field<T>::N() const {
    return _N;
}



template<class T>
size_t Field<T>::N(const size_t i) const {
    return _N[i];
}



template<class T>
std::array<double, 3> Field<T>::L() const {
    return _L;
}



template<class T>
double Field<T>::L(const size_t i) const {
    return _L[i];
}
