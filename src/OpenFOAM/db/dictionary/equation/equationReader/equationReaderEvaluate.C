/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::equationReader::evaluate
(
    const word& equationName
)
{
    label index(lookup(equationName));
    if (index < 0)
    {
        FatalErrorIn("equationReader::evaluate(const word)")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    return evaluate(index);
}


Foam::dimensionedScalar Foam::equationReader::evaluate
(
    const label& index,
    label storageOffset
)
{
    if (debug)
    {
        Info << "Evaluating equation " << eqns_[index].equationName()
            << " at index " << index << ", given by:" << token::NL
            << token::TAB << eqns_[index].rawText() << endl;
    }
    if ((index < 0) || (index >= eqns_.size()))
    {
        FatalErrorIn("equationReader::update(const index)")
            << "Index " << index << " out of bounds (0, " << eqns_.size() - 1
            << ")"
            << abort(FatalError);
    }
    
    if (eqns_[index].size() == 0)
    {
        parse(index);
    }

    label storeIndex(-1);
    dimensionedScalar ds("empty", dimless, 0);
    
    // Bug fix - dropped reference
    // equationOperationList& ops(this->operator[](index));

    for (label i(0); i < eqns_[index].size(); i++)
    {
        if
        (
            (ds.name() == "empty")
         && (
                eqns_[index].ops()[i].operation()
             != equationOperation::otretrieve
            )
        )
        {
            FatalErrorIn("equationReader::update(index)")
                << "Bad operation list.  Operation at " << i << " either "
                << "follows a 'store', or is the first operation.  Therefore "
                << "it should be retrieve, but it is "
                << eqns_[index].ops()[i].operation() << "."
                << abort(FatalError);
        }
        
        dimensionedScalar source("noSource", dimless, 0);
        
        // Special case - only otstore does not require source
        if
        (
            eqns_[index].ops()[i].operation()
         != equationOperation::otstore
        )
        {
            dsEqual
            (
                source,
                getSource
                (
                    index,
                    i,
                    storeIndex + storageOffset,
                    storageOffset
                )
            );
        }

        if (debug > 1)
        {
            Info << "Performing operation: ["
                << equationOperation::opName
                (
                    eqns_[index].ops()[i].operation()
                ) << "] using source [";
            if
            (
                eqns_[index].ops()[i].sourceList()
             == equationOperation::sldictSource
            )
            {
                Info << dictLookups_[eqns_[index].ops()[i].dictLookupIndex()];
            }
            else if
            (
                eqns_[index].ops()[i].sourceList()
             == equationOperation::slequation
            )
            {
                Info << eqns_
                    [
                        mag(eqns_[index].ops()[i].sourceIndex() - 1)
                    ].equationName();
            }
            else if
            (
                eqns_[index].ops()[i].sourceList()
             == equationOperation::slexternalScalar
            )
            {
                Info << externalScalarNames_
                    [
                        mag(eqns_[index].ops()[i].sourceIndex() - 1)
                    ];
            }
            else
            {
                Info << ds.name();
            }
            Info << "] read from ["
                << equationOperation::sourceName
                (
                    eqns_[index].ops()[i].sourceList()
                ) << "]..." << endl;
        }
        switch (eqns_[index].ops()[i].operation())
        {
            case equationOperation::otnone:
                FatalErrorIn("equationReader::update(index)")
                    << "Empty operation called.  Empty operations should only "
                    << "exist temporarily during parsing, and they should not "
                    << "remain in the operation list at this point.  Either "
                    << "you have corrupt data, or this is a bug."
                    << abort(FatalError);
                break;
            case equationOperation::otretrieve:
                dsEqual(ds, source);
                break;
            case equationOperation::otstore:
                storeIndex++;
                storage_.setSize(storeIndex + storageOffset + 1);
                storage_.set
                (
                    storeIndex + storageOffset,
                    new dimensionedScalar(ds)
                );
                dsEqual(ds, dimensionedScalar("empty", dimless, 0));
                break;
            case equationOperation::otplus:
                if (eqns_[index].changeDimensions())
                {
                    dimensionedScalar dsTemp
                    (
                        "dsTemp",
                        dimless,
                        ds.value() + source.value()
                    );
                    dsEqual(ds, dsTemp);
                }
                else
                {
                    if
                    (
                        dimensionSet::debug
                     && (ds.dimensions() != source.dimensions())
                    )
                    {
                        WarningIn("equationReader::evaluate")
                            << "Dimension error thrown for equation "
                            << eqns_[index].equationName() << ", given by:"
                            << token::NL << token::TAB
                            << eqns_[index].rawText();
                    }
                    ds = ds + source;
                }
                break;
            case equationOperation::otminus:
                if (eqns_[index].changeDimensions())
                {
                    dimensionedScalar dsTemp
                    (
                        "dsTemp",
                        dimless,
                        ds.value() - source.value()
                    );
                    dsEqual(ds, dsTemp);
                }
                else
                {
                    if
                    (
                        dimensionSet::debug
                     && (ds.dimensions() != source.dimensions())
                    )
                    {
                        WarningIn("equationReader::evaluate")
                            << "Dimension error thrown for equation "
                            << eqns_[index].equationName() << ", given by:"
                            << token::NL << token::TAB
                            << eqns_[index].rawText() << token::NL << endl;
                    }
                    ds = ds - source;
                }
                break;
            case equationOperation::ottimes:
                dsEqual(ds, ds * source);
                break;
            case equationOperation::otdivide:
                dsEqual(ds, ds / source);
                break;
            case equationOperation::otpow:
                dsEqual(ds, pow(ds, source));
                break;
            case equationOperation::otsign:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'sign' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                ds = sign(ds);
                break;
            case equationOperation::otpos:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'pos' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                ds = pos(ds);
                break;
            case equationOperation::otneg:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'neg' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                ds = neg(ds);
                break;
            case equationOperation::otmag:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'mag' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                ds = mag(ds);
                break;
            case equationOperation::otlimit:
                if (source.name() == "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'limit' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                ds.value() = limit(ds.value(), source.value());
                break;
            case equationOperation::otminMod:
                if (source.name() == "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'minMod' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                ds.value() = minMod(ds.value(), source.value());
                break;
            case equationOperation::otsqrtSumSqr:
                if (source.name() == "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'sqrtSumSqr' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                ds.value() = sqrtSumSqr(ds.value(), source.value());
                break;
            case equationOperation::otsqr:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'sqr' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, sqr(ds));
                break;
            case equationOperation::otpow3:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'pow3' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, pow3(ds));
                break;
            case equationOperation::otpow4:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'pow4' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, pow4(ds));
                break;
            case equationOperation::otpow5:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'pow5' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, pow5(ds));
                break;
            case equationOperation::otinv:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'inv' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, pTraits<scalar>::one / ds);
                break;
            case equationOperation::otsqrt:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'sqrt' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, sqrt(ds));
                break;
            case equationOperation::otcbrt:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'cbrt' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, cbrt(ds));
                break;
            case equationOperation::othypot:
                if (source.name() == "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'hypot' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                dsEqual(ds, hypot(ds, source));
            case equationOperation::otexp:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'exp' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, exp(ds));
                break;
            case equationOperation::otlog:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'log' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, log(ds));
                break;
            case equationOperation::otlog10:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'log10' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, log10(ds));
                break;
            case equationOperation::otsin:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'sin' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, sin(ds));
                break;
            case equationOperation::otcos:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'cos' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, cos(ds));
                break;
            case equationOperation::ottan:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'tan' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, tan(ds));
                break;
            case equationOperation::otasin:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'asin' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, asin(ds));
                break;
            case equationOperation::otacos:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'acos' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, acos(ds));
                break;
            case equationOperation::otatan:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'atan' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, atan(ds));
                break;
            case equationOperation::otsinh:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'sinh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, sinh(ds));
                break;
            case equationOperation::otcosh:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'cosh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, cosh(ds));
                break;
            case equationOperation::ottanh:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'tanh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, tanh(ds));
                break;
            case equationOperation::otasinh:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'asinh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, sinh(ds));
                break;
            case equationOperation::otacosh:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'acosh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, cosh(ds));
                break;
            case equationOperation::otatanh:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'atanh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, tanh(ds));
                break;
            case equationOperation::oterf:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'erf' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, erf(ds));
                break;
            case equationOperation::oterfc:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'erfc' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, erfc(ds));
                break;
            case equationOperation::otlgamma:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'lgamma' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, lgamma(ds));
                break;
            case equationOperation::otj0:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'j0' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, j0(ds));
                break;
            case equationOperation::otj1:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'j1' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, j1(ds));
                break;
            case equationOperation::otjn:
            {
                if (source.name() == "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'jn' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                int di(ds.value());
                dsEqual(ds, jn(di, source));
                break;
            }
            case equationOperation::oty0:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'y0' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, y0(ds));
                break;
            case equationOperation::oty1:
                if (source.name() != "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'y1' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                dsEqual(ds, y1(ds));
                break;
            case equationOperation::otyn:
            {
                if (source.name() == "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'yn' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                int di(ds.value());
                dsEqual(ds, yn(di, source));
                break;
            }
            case equationOperation::otmax:
            {
                if (source.name() == "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'max' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                dsEqual(ds, max(ds, source));
                break;
            }
            case equationOperation::otmin:
            {
                if (source.name() == "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'min' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                dsEqual(ds, min(ds, source));
                break;
            }
            case equationOperation::otstabilise:
            {
                if (source.name() == "noSource")
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqns_[index].equationName() << ", given by "
                        << token::NL << token::TAB  << eqns_[index].rawText()
                        << token::NL << "Function 'stabilise' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                ds.value() = stabilise(ds.value(), source.value());
                break;
            } // end case
        } // end switch
        if (debug > 1)
        {
            Info << "Operaion result is " << ds << endl;
        }
    }

    ds.name() = eqns_[index].equationName();
    
    //Move one level back up on the dependents_ list
    if (dependents_.size())
    {
        dependents_.setSize(dependents_.size() - 1);
    }    

    storage_.setSize(storageOffset);
    if (eqns_[index].changeDimensions())
    {
        ds.dimensions().reset(eqns_[index].overrideDimensions());
    }
    if (debug)
    {
        Info << "Equation evaluated.  Result is: " << ds << endl;
    }
    dsEqual(eqns_[index].lastResult(), ds);
    return ds;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
