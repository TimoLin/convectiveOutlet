/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "convectiveOutletFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    convecSpeed_(0.0),
    fieldInf_(pTraits<Type>::zero),
    lInf_(-GREAT)
{
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const convectiveOutletFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    convecSpeed_(ptf.convecSpeed_),
    fieldInf_(ptf.fieldInf_),
    lInf_(ptf.lInf_)
{}


template<class Type>
Foam::convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    convecSpeed_(dict.lookupOrDefault<scalar>("convecSpeed",0.0)),
    fieldInf_(pTraits<Type>::zero),
    lInf_(-GREAT)
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    this->refValue() = *this;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;

    if (dict.readIfPresent("lInf", lInf_))
    {
        dict.lookup("fieldInf") >> fieldInf_;

        if (lInf_ < 0.0)
        {
            FatalIOErrorIn
            (
                "convectiveOutletFvPatchField<Type>::"
                "convectiveOutletFvPatchField"
                "("
                    "const fvPatch&, "
                    "const DimensionedField<Type, volMesh>&, "
                    "const dictionary&"
                ")",
                dict
            )   << "unphysical lInf specified (lInf < 0)" << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << exit(FatalIOError);
        }
    }
}


template<class Type>
Foam::convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const convectiveOutletFvPatchField& ptpsf
)
:
    mixedFvPatchField<Type>(ptpsf),
    phiName_(ptpsf.phiName_),
    rhoName_(ptpsf.rhoName_),
    fieldInf_(ptpsf.fieldInf_),
    lInf_(ptpsf.lInf_)
{}


template<class Type>
Foam::convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const convectiveOutletFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptpsf, iF),
    phiName_(ptpsf.phiName_),
    rhoName_(ptpsf.rhoName_),
    convecSpeed_(ptpsf.convecSpeed_),
    fieldInf_(ptpsf.fieldInf_),
    lInf_(ptpsf.lInf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::convectiveOutletFvPatchField<Type>::convectionSpeed() const
{
    if (convecSpeed_ == 0.0)
    {
        const surfaceScalarField& phi =
            this->db().objectRegistry::template lookupObject<surfaceScalarField>
            (phiName_);

        fvsPatchField<scalar> phip =
            this->patch().template lookupPatchField<surfaceScalarField, scalar>
            (
                phiName_
            );

        if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
        {
            // Compressible case
            const fvPatchScalarField& rhop =
                this->patch().template lookupPatchField<volScalarField, scalar>
                (
                    rhoName_
                );

            tmp<scalarField> Un(
                    new scalarField
                    (
                        phip.size(),
                        mag(gSum(phip)/gSum(rhop*this->patch().magSf()))                        
                    )
                    );
            return Un;
        }
        else
        {
            // Incompressible case
            // Calculate outlet patch mean velocity
            tmp<scalarField> Un(
                    new scalarField
                    (
                        phip.size(),
                        mag(gSum(phip)/gSum(this->patch().magSf()))                        
                    )
                    );
            return Un;
        }
    }
    else
    {
        fvsPatchField<scalar> phip =
            this->patch().template lookupPatchField<surfaceScalarField, scalar>
            (
             phiName_
            );
        tmp<scalarField> Un(new scalarField(phip.size(),convecSpeed_));
        return Un;
    }
}


template<class Type>
void Foam::convectiveOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    word ddtScheme
    (
        this->dimensionedInternalField().mesh()
       .ddtScheme(this->dimensionedInternalField().name())
    );
    scalar deltaT = this->db().time().deltaTValue();

    const GeometricField<Type, fvPatchField, volMesh>& field =
        this->db().objectRegistry::template
        lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            this->dimensionedInternalField().name()
        );

    // Calculate the convection speed of the field wave
    // If the wave is incoming set the speed to 0.
    scalarField w(convectionSpeed());//, scalar(0)));

    // Output convectionSpeed to stdout
    Info<<"Outlet bondary convection speed: "<<w[0]<<" m/s"<<endl;

    // Calculate the field wave coefficient alpha (See notes)
    //  deltaCoeffs() is the reciprocal distance
    //    Between owner and neighbour cell centroids of an internal face
    const scalarField alpha(w*deltaT*this->patch().deltaCoeffs());

    label patchi = this->patch().index();

    // Non-reflecting outflow boundary
    // If lInf_ defined setup relaxation to the value fieldInf_.
    if (lInf_ > 0)
    {
        // Calculate the field relaxation coefficient k (See notes)
        const scalarField k(w*deltaT/lInf_);

        if
        (
            ddtScheme == fv::EulerDdtScheme<scalar>::typeName
         || ddtScheme == fv::CrankNicolsonDdtScheme<scalar>::typeName
        )
        {
            this->refValue() =
            (
                field.oldTime().boundaryField()[patchi] + k*fieldInf_
            )/(1.0 + k);

            this->valueFraction() = (1.0 + k)/(1.0 + alpha + k);
        }
        else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
        {
            this->refValue() =
            (
                2.0*field.oldTime().boundaryField()[patchi]
              - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
              + k*fieldInf_
            )/(1.5 + k);

            this->valueFraction() = (1.5 + k)/(1.5 + alpha + k);
        }
        else
        {
            FatalErrorIn("convectiveOutletFvPatchField<Type>::updateCoeffs()")
                << "    Unsupported temporal differencing scheme : "
                << ddtScheme << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            ddtScheme == fv::EulerDdtScheme<scalar>::typeName
         || ddtScheme == fv::CrankNicolsonDdtScheme<scalar>::typeName
        )
        {
            this->refValue() = field.oldTime().boundaryField()[patchi];

            this->valueFraction() = 1.0/(1.0 + alpha);
        }
        else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
        {
            this->refValue() =
            (
                2.0*field.oldTime().boundaryField()[patchi]
              - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
            )/1.5;

            this->valueFraction() = 1.5/(1.5 + alpha);
        }
        else
        {
            FatalErrorIn
            (
                "convectiveOutletFvPatchField<Type>::updateCoeffs()"
            )   << "    Unsupported temporal differencing scheme : "
                << ddtScheme
                << "\n    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << exit(FatalError);
        }
    }
    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::convectiveOutletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    this->template writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    this->template writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

    if (lInf_ > 0)
    {
        os.writeKeyword("fieldInf") << fieldInf_ << token::END_STATEMENT << nl;
        os.writeKeyword("lInf") << lInf_ << token::END_STATEMENT << nl;
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
