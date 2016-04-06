// --------------------- //
//    Static Members     //
// --------------------  //

template <typename Base, const char* file_suffix>
unsigned int Tracer<Base, file_suffix>::size = sizeof(typename Base::RealType);

template <typename Base, const char* file_suffix>
unsigned int Tracer<Base, file_suffix>::total_size = 0;

template <typename Base, const char* file_suffix>
unsigned int Tracer<Base, file_suffix>::object_count = 0;

template <typename Base, const char* file_suffix>
std::vector<unsigned int> Tracer<Base, file_suffix>::object_counts = std::vector<unsigned int>();

template <typename Base, const char* file_suffix>
std::vector<unsigned int> Tracer<Base, file_suffix>::total_sizes = std::vector<unsigned int>();

template <typename Base, const char* file_suffix>
void Tracer<Base, file_suffix>::start_tracing()
{
    object_count  = 0;
    total_size    = 0;
    object_counts = std::vector<unsigned int>();
    total_sizes   = std::vector<unsigned int>();
}

template <typename Base, const char* file_suffix>
void Tracer<Base, file_suffix>::stop_tracing(const std::string &prefix)
{
    std::string filename = prefix;
    filename += "_";
    filename += file_suffix;
    filename += ".dat";

    std::ofstream file;
    file.open(filename);

    for (unsigned int i = 0; i < total_sizes.size(); i++)
    {
        file << object_counts[i] << " " << total_sizes[i] << std::endl;
    }

    file.close();
}

// --------------------- //
//    Member Functions   //
// --------------------  //

template <typename Base, const char* file_suffix>
Tracer<Base, file_suffix>::Tracer()
    : Base()
{
    object_count++;
    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base, const char* file_suffix>
Tracer<Base, file_suffix>::Tracer(const Tracer &t)
    : Base(t)
{
    object_count++;
    total_size += this->cols() * this->rows() * size;

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base, const char* file_suffix>
Tracer<Base, file_suffix>::Tracer(const Base &t)
    : Base(t)
{
    object_count++;
    total_size += this->cols() * this->rows() * size;

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base, const char* file_suffix>
Tracer<Base, file_suffix>::Tracer(Tracer &&t)
    : Base(std::move(t))
{
    object_count++;

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base, const char* file_suffix>
Tracer<Base, file_suffix>::Tracer(Base &&t)
    : Base(std::move(t))
{
    object_count++;
    total_size += this->cols() * this->rows() * size;

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base, const char* file_suffix>
Tracer<Base, file_suffix> & Tracer<Base, file_suffix>::operator=(const Tracer &t)
{
    total_size -= this->cols() * this->rows() * size;
    (*this).Base::operator=(t);
    total_size += this->cols() * this->rows() * size;

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);

    return *this;
}

template <typename Base, const char* file_suffix>
Tracer<Base, file_suffix> & Tracer<Base, file_suffix>::operator=(const Base &t)
{
    total_size -= this->cols() * this->rows() * size;
    (*this).Base::operator=(t);
    total_size += this->cols() * this->rows() * size;

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);

    return *this;
}

template <typename Base, const char* file_suffix>
Tracer<Base, file_suffix> & Tracer<Base, file_suffix>::operator=(Tracer &&t)
{
    total_size -= this->cols() * this->rows() * size;
    (*this).Base::operator=(std::move(t));

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);

    return *this;
}

template <typename Base, const char* file_suffix>
Tracer<Base, file_suffix> & Tracer<Base, file_suffix>::operator=(Base &&t)
{
    total_size -= this->cols() * this->rows() * size;
    (*this).Base::operator=(std::move(t));
    total_size += this->cols() * this->rows() * size;

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);

    return *this;
}

template <typename Base, const char* file_suffix>
Tracer<Base, file_suffix>::~Tracer()
{
    total_size -= this->cols() * this->rows() * size;
    object_count--;

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base, const char* file_suffix>
    template<typename T1>
void Tracer<Base, file_suffix>::resize(T1 m, T1 n)
{
    total_size -= this->cols() * this->rows() * size;
    total_size += m * n * size;
    (*this).Base::resize(m, n);

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base, const char* file_suffix>
    template<typename T1>
void Tracer<Base, file_suffix>::resize(T1 m)
{
    total_size -= this->cols() * size;
    total_size += m * size;
    (*this).Base::resize(m);

    object_counts.push_back(object_count);
    total_sizes.push_back(total_size);
}

template <typename Base, const char* file_suffix>
    template<typename T1>
unsigned int Tracer<Base, file_suffix>::cols_base(decltype(&T1::cols)) const
{
    return this->Base::cols();
}

template <typename Base, const char* file_suffix>
    template<typename T1>
unsigned int Tracer<Base, file_suffix>::cols_base(...) const
{
    return 1;
}

template <typename Base, const char* file_suffix>
unsigned int Tracer<Base, file_suffix>::cols() const
{
    return cols_base<Base>(0);
}
