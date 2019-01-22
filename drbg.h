#ifndef DRBG_H
#define DRBG_H

#include <QSharedPointer>
#include <QString>
#include <QCryptographicHash>

class DRBG
{
    struct pimpl;
    QSharedPointer<pimpl> _pimpl;

public:
    DRBG(const QByteArray &seed_material = QByteArray(), int entropy = 0);

    static DRBG fromRandomSeed();

    int entropy() const;
    QByteArray seed() const;

    QByteArray generate(long cb) const;
    void reset();
};

#endif // DRBG_H
