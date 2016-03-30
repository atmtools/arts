// --------------------- //
//    Static Members     //
// --------------------  //

template <typename Base>
unsigned int Tracer<Base>::size = sizeof(typename Base::RealType);

template <typename Base>
unsigned int Tracer<Base>::total_size = 0;

template <typename Base>
unsigned int Tracer<Base>::object_count = 0;

template <typename Base>
std::vector<unsigned int> Tracer<Base>::object_counts = std::vector<unsigned int>();

template <typename Base>
std::vector<unsigned int> Tracer<Base>::total_sizes = std::vector<unsigned int>();

template <typename Base>
void Tracer<Base>::start_tracing()
{
    object_count  = 0;
    total_size    = 0;
    object_counts = std::vector<unsigned int>();
    total_sizes   = std::vector<unsigned int>();
}

template <typename Base>
void Tracer<Base>::stop_tracing(const std::string &prefix)
{
    std::string filename = prefix;
    filename += "_";
    filename += TypeName<Base>::name;
    filename += ".dat";

    std::ofstream file;
    file.open(filename);

    for (unsigned int i = 0; i < object_counts.size(); i++)
    {
        file << object_counts[i] << " " << total_sizes[i] << std::endl;
    }

    file.close();
}

// --------------------- //
//    Member Functions   //
// --------------------  //

template <typename Base>
    template <typename T>
Tracer<Base>::Tracer(T &&t)
    : Base(std::forward<T>(t))
{
    object_count++;
    total_size += this->cols() * this->rows() * size;

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base>
    template <typename T>
Tracer<Base> & Tracer<Base>::operator=(T &&t)
{
    total_size -= this->cols() * this->rows() * size;
    (*this).Base::operator=(std::forward<T>(t));
    total_size += this->cols() * this->rows() * size;

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base>
Tracer<Base>::~Tracer()
{
    total_size -= this->cols() * this->rows() * size;
    object_count--;

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base>
    template<typename T1>
void Tracer<Base>::resize(T1 m, T1 n)
{
    total_size -= this->cols() * this->rows() * size;
    total_size += m * n * size;
    (*this).Base::resize(m, n);

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base>
    template<typename T1>
void Tracer<Base>::resize(T1 m)
{
    total_size -= this->cols() * size;
    total_size += m * size;
    (*this).Base::resize(m);

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base>
    template<typename T1>
unsigned int Tracer<Base>::cols_base(decltype(&T1::cols)) const
{
    return this->Base::cols();
}

template <typename Base>
    template<typename T1>
unsigned int Tracer<Base>::cols_base(...) const
{
    return 1;
}

template <typename Base>
unsigned int Tracer<Base>::cols() const
{
    return cols_base<Base>(0);
}
