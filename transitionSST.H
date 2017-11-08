/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

Class
    Foam::RASModels::transitionSST

Group
    grpRASTurbulence

See also
    Foam::transitionSST

SourceFiles
    transitionSST.C

\*---------------------------------------------------------------------------*/

#ifndef transitionSST_H
#define transitionSST_H

// #include "kOmegaSSTBase.H"
#include "kOmegaSST.H"
// #include "RASModel.H"
// #include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class transitionSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class transitionSST
:
    public kOmegaSST<BasicTurbulenceModel>
{
    // Disallow default bitwise copy construct and assignment
        transitionSST(const transitionSST&);
        void operator=(const transitionSST&);

protected:
    // volScalarField k_;
    // volScalarField omega_;
    volScalarField ReThetatTilda_;
    volScalarField gamma_;
    // volScalarField::Internal gammaEff_;

    tmp<volScalarField> Pgamma() const;
    tmp<volScalarField> Fonset1() const;
    tmp<volScalarField> Fonset2() const;
    tmp<volScalarField> Fonset3() const;
    tmp<volScalarField> Fonset() const;
    tmp<volScalarField> ReV() const;
    tmp<volScalarField> Rt() const;
    tmp<volScalarField> Flength() const;
    tmp<volScalarField> ReThetac() const;
    tmp<volScalarField> Fsublayer() const;
    tmp<volScalarField> ReOmega() const;
    tmp<volScalarField> gammaSep() const;
    tmp<volScalarField> Freattach() const;
    tmp<volScalarField> gammaEff() const;
    tmp<volScalarField> Fturb() const;

    tmp<volScalarField> PReThetatTilda() const;
    tmp<volScalarField> Fthetat() const;
    tmp<volScalarField> ReThetat() const;
    scalar F(scalar lambda, scalar Tu, scalar dUds) const;
    tmp<volScalarField> delta() const;
    tmp<volScalarField> Fwake() const;

    tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
    tmp<volScalarField> F3() const;
    tmp<volScalarField> Ry() const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("transitionSST");

    // Constructors

        //- Construct from components
        transitionSST
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~transitionSST()
    {}

    tmp<volScalarField> DgammaEff() const
    {
        return this->mu() + this->mut();
    }

    tmp<volScalarField> DReThetatTildaEff() const
    {
        return 2.0*(this->mu() + this->mut());
    }

    virtual tmp<volScalarField> Pk(const volScalarField G) const;
    virtual tmp<volScalarField> epsilonByk (const volScalarField & F1, const volScalarField& F2) const;
    virtual void correctGammaReThetatTilda();
    
    const volScalarField& gamma() const
    {
        return gamma_;
    }

    const volScalarField& ReThetatTilda() const
    {
        return ReThetatTilda_;
    }
    
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "transitionSST.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
