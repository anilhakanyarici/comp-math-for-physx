#include "drbg.h"

struct DRBG::pimpl
{
    mutable int _position;
    int _entropy;
    QCryptographicHash::Algorithm _hash_alg;
    QByteArray _seed_material;
    QByteArray _current_seed;
    mutable QByteArray _last_hash;
};

DRBG::DRBG(const QByteArray &seed_material, int entropy)
{
    assert(seed_material.size() < 256 && "Length of seedMaterial cannot be longer than 255.");
    assert(entropy >= 0 && "Entropy cannot be smaller than zero.");

    this->_pimpl = QSharedPointer<pimpl>(new pimpl);
    this->_pimpl->_hash_alg = QCryptographicHash::Sha512;
    this->_pimpl->_seed_material = seed_material;
    QCryptographicHash hash(this->_pimpl->_hash_alg);
    hash.addData(QByteArray().append(::clock()));
    this->_pimpl->_entropy = entropy;
    QByteArray entropy_bytes = hash.result().left(this->_pimpl->_entropy);
    this->_pimpl->_current_seed = entropy_bytes.append(seed_material);
    hash.reset();
    hash.addData(this->_pimpl->_current_seed);
    this->_pimpl->_last_hash = hash.result();
    this->_pimpl->_position = 0;
}

DRBG DRBG::fromRandomSeed()
{
    QCryptographicHash hash(QCryptographicHash::Sha512);
    hash.addData(QByteArray().append(::clock()));
    return DRBG(hash.result());
}

int DRBG::entropy() const
{
    return this->_pimpl->_entropy;
}

QByteArray DRBG::seed() const
{
    return this->_pimpl->_seed_material;
}

QByteArray DRBG::generate(long cb) const
{
    if (cb == 0)
        return QByteArray();
    assert(cb >= 0 && "cb cannot be smaller than 0.");

    int readable = 64 - this->_pimpl->_position;

    if (readable >= cb)
    {
        QByteArray cbBytes(cb, Qt::Initialization::Uninitialized);
        ::memcpy(cbBytes.data(), this->_pimpl->_last_hash.data() + this->_pimpl->_position, cb);
        this->_pimpl->_position += cb;
        return cbBytes;
    }
    else
    {
        QCryptographicHash hash(this->_pimpl->_hash_alg);
        if (readable == 0)
        {
            hash.addData(QByteArray().append(::clock()));
            QByteArray entropy_bytes = hash.result().left(this->_pimpl->_entropy);
            hash.reset();
            hash.addData(this->_pimpl->_last_hash.append(entropy_bytes));
            this->_pimpl->_last_hash = hash.result();
            readable = 64;
            this->_pimpl->_position = 0;
        }
        QByteArray cbBytes(cb, Qt::Initialization::Uninitialized);
        if (readable >= cb)
        {
            ::memcpy(cbBytes.data(), this->_pimpl->_last_hash.data() + this->_pimpl->_position, cb);
            this->_pimpl->_position += cb;
            return cbBytes;
        }
        ::memcpy(cbBytes.data(), this->_pimpl->_last_hash.data() + this->_pimpl->_position, cb);
        cb -= readable;
        long cbOffset = readable;
        this->_pimpl->_position = 0;

        hash.reset();
        hash.addData(QByteArray().append(::clock()));
        QByteArray entropy_bytes = hash.result().left(this->_pimpl->_entropy);
        hash.reset();
        hash.addData(this->_pimpl->_last_hash.append(entropy_bytes));
        this->_pimpl->_last_hash = hash.result();

        while (cb > 63)
        {
            ::memcpy(cbBytes.data() + cbOffset, this->_pimpl->_last_hash.data(), 64);
            cb -= 64;
            cbOffset += 64;

            hash.reset();
            hash.addData(QByteArray().append(::clock()));
            QByteArray entropy_bytes = hash.result().left(this->_pimpl->_entropy);
            hash.reset();
            hash.addData(this->_pimpl->_last_hash.append(entropy_bytes));
            this->_pimpl->_last_hash = hash.result();
        }
        if (cb > 0)
        {
            ::memcpy(cbBytes.data() + cbOffset, this->_pimpl->_last_hash.data(), cb);
            this->_pimpl->_position += cb;
        }
        return cbBytes;
    }
}

void DRBG::reset()
{
    this->_pimpl->_hash_alg = QCryptographicHash::Sha512;
    QByteArray seed_material = this->_pimpl->_seed_material;
    QCryptographicHash hash(this->_pimpl->_hash_alg);
    hash.addData(QByteArray().append(::clock()));
    QByteArray entropy_bytes = hash.result().left(this->_pimpl->_entropy);
    this->_pimpl->_current_seed = entropy_bytes.append(seed_material);
    hash.reset();
    hash.addData(this->_pimpl->_current_seed);
    this->_pimpl->_last_hash = hash.result();
    this->_pimpl->_position = 0;
}
