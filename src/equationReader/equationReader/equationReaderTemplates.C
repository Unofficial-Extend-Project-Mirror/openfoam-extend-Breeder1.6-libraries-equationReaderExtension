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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::equationReader::evaluateTypeField
(
    Field<Type>& resultField,
    const word& componentName,
    const word& equationName,
    const label geoIndex
) const
{
    // Get equationIndex
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateTypeField")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }

    // Get componentIndex
    Type dummy;
    label componentIndex(-1);
    forAll(dummy, i)
    {
        if (Type::componentNames[i] == componentName)
        {
            componentIndex = i;
            break;
        }
    }
    if (componentIndex < 0)
    {
        wordList validNames(dummy.size());
        forAll(dummy, i)
        {
            validNames[i] = Type::componentNames[i];
        }
        FatalErrorIn("equationReader::evaluateTypeField")
            << componentName << " is not a valid component name.  Valid names "
            << "are " << validNames
            << abort(FatalError);
    }

    evaluateTypeField(resultField, componentIndex, equationIndex, geoIndex);
}


template<class Type>
void Foam::equationReader::evaluateTypeField
(
    Field<Type>& resultField,
    const label componentIndex,
    const label equationIndex,
    const label geoIndex
) const
{
    forAll(resultField, cellIndex)
    {
        resultField[cellIndex][componentIndex] =
            evaluateScalar(equationIndex, cellIndex, geoIndex);
    }
}


template<class GeoMesh>
void Foam::equationReader::evaluateDimensionedScalarField
(
    DimensionedField<scalar, GeoMesh>& resultDField,
    const word& equationName,
    const label geoIndex
) const
{
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateDimensionedScalarField")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    evaluateDimensionedScalarField(resultDField, equationIndex, geoIndex);
}


template<class GeoMesh>
void Foam::equationReader::evaluateDimensionedScalarField
(
    DimensionedField<scalar, GeoMesh>& resultDField,
    const label equationIndex,
    const label geoIndex
) const
{
    forAll(resultDField.field(), cellIndex)
    {
        resultDField.field()[cellIndex] =
            evaluateScalar(equationIndex, cellIndex, geoIndex);
        checkFinalDimensions
        (
            equationIndex,
            resultDField.dimensions(),
            resultDField.name()
        );
    }
}


template<class Type, class GeoMesh>
void Foam::equationReader::evaluateDimensionedTypeField
(
    DimensionedField<Type, GeoMesh>& resultDField,
    const word& componentName,
    const word& equationName,
    const label geoIndex
) const
{
    // Get equationIndex
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateDimensionedTypeField")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }

    // Get componentIndex
    Type dummy;
    label componentIndex(-1);
    forAll(dummy, i)
    {
        if (Type::componentNames[i] == componentName)
        {
            componentIndex = i;
            break;
        }
    }
    if (componentIndex < 0)
    {
        wordList validNames(dummy.size());
        forAll(dummy, i)
        {
            validNames[i] = Type::componentNames[i];
        }
        FatalErrorIn("equationReader::evaluateDimensionedTypeField")
            << componentName << " is not a valid component name.  Valid names "
            << "are " << validNames
            << abort(FatalError);
    }

    evaluateDimensionedTypeField
    (
        resultDField,
        componentIndex,
        equationIndex,
        geoIndex
    );
}


template<class Type, class GeoMesh>
void Foam::equationReader::evaluateDimensionedTypeField
(
    DimensionedField<Type, GeoMesh>& resultDField,
    const label componentIndex,
    const label equationIndex,
    const label geoIndex
) const
{
    forAll(resultDField.field(), cellIndex)
    {
        resultDField.field()[cellIndex][componentIndex] =
            evaluateScalar(equationIndex, cellIndex, geoIndex);
        checkFinalDimensions
        (
            equationIndex,
            resultDField.dimensions(),
            resultDField.name() + "." + Type::componentNames[componentIndex]
        );
    }
}


template<template<class> class PatchField, class GeoMesh>
void Foam::equationReader::evaluateGeometricScalarField
(
    GeometricField<scalar, PatchField, GeoMesh>& resultGField,
    const word& equationName
) const
{
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateGeometricScalarField")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    evaluateGeometricScalarField(resultGField, equationIndex);
}


template <template<class> class PatchField, class GeoMesh>
void Foam::equationReader::evaluateGeometricScalarField
(
    GeometricField<scalar, PatchField, GeoMesh>& resultGField,
    const label equationIndex
) const
{
    // Internal field: geoIndex = 0
    forAll(resultGField.internalField(), cellIndex)
    {
        resultGField.internalField()[cellIndex] =
            evaluateScalar(equationIndex, cellIndex, 0);
    }
    // Boundary fields: geoIndex = patchIndex + 1
    forAll(resultGField.boundaryField(), patchIndex)
    {
        forAll(resultGField.boundaryField()[patchIndex], cellIndex)
        {
            resultGField.boundaryField()[patchIndex][cellIndex] =
                evaluateScalar
                (
                    equationIndex,
                    cellIndex,
                    patchIndex + 1
                );
        }
    }
    checkFinalDimensions
    (
        equationIndex,
        resultGField.dimensions(),
        resultGField.name()
    );
}


template
<
    class Type, template<class> class PatchField, class GeoMesh
>
void Foam::equationReader::evaluateGeometricTypeField
(
    GeometricField<Type, PatchField, GeoMesh>& resultGField,
    const word& componentName,
    const word& equationName
) const
{
    // Get equationIndex
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateGeometricTypeField")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }

    // Get componentIndex
    Type dummy;
    label componentIndex(-1);
    forAll(dummy, i)
    {
        if (Type::componentNames[i] == componentName)
        {
            componentIndex = i;
            break;
        }
    }
    if (componentIndex < 0)
    {
        wordList validNames(dummy.size());
        forAll(dummy, i)
        {
            validNames[i] = Type::componentNames[i];
        }
        FatalErrorIn("equationReader::evaluateGeometricTypeField")
            << componentName << " is not a valid component name.  Valid names "
            << "are " << validNames
            << abort(FatalError);
    }
    evaluateGeometricTypeField(resultGField, componentIndex, equationIndex);
}


template
<
    class Type, template<class> class PatchField, class GeoMesh
>
void Foam::equationReader::evaluateGeometricTypeField
(
    GeometricField<Type, PatchField, GeoMesh>& resultGField,
    const label componentIndex,
    const label equationIndex
) const
{
    // Internal field: geoIndex = 0
    forAll(resultGField.internalField(), cellIndex)
    {
        resultGField.internalField()[cellIndex][componentIndex] =
            evaluateScalar(equationIndex, cellIndex, 0);
    }
    // Boundary fields: geoIndex = patchIndex + 1
    forAll(resultGField.boundaryField(), patchIndex)
    {
        forAll(resultGField.boundaryField()[patchIndex], cellIndex)
        {
            resultGField.boundaryField()
                [patchIndex]
                [cellIndex]
                [componentIndex]
              = evaluateScalar
                (
                    equationIndex,
                    cellIndex,
                    patchIndex + 1
                );
        }
    }
    checkFinalDimensions
    (
        equationIndex,
        resultGField.dimensions(),
        resultGField.name() + "." + Type::componentNames[componentIndex]
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

