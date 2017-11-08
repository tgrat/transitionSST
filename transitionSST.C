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

\*---------------------------------------------------------------------------*/

#include "transitionSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
transitionSST<BasicTurbulenceModel>::transitionSST
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSST<BasicTurbulenceModel>
    (
        // type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    ReThetatTilda_
    (
        IOobject
        (
            IOobject::groupName("ReThetatTilda", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    gamma_
    (
        IOobject
        (
            IOobject::groupName("gamma", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )

    //  gammaEff_
    // (
    //     IOobject
    //     (
    //         IOobject::groupName("gammaEff", U.group()),
    //         this->runTime_.timeName(),
    //         this->mesh_,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     this->mesh_
    // )

{
    // if (type == typeName)
    // {
    //     this->printCoeffs(type);
    // }
}

// * * * * * * * * * * * * * * Altered member functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Pgamma() const
{
    return 2.0*Flength()*this->rho_*sqrt(2.0*magSqr(symm(fvc::grad(this->U()))))*sqrt(gamma_*Fonset())*(1.0 - gamma_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Fonset1() const
{
    return ReV()/(scalar(2.193)*ReThetac());
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Fonset2() const
{
    return min(
                max(
                    Fonset1(),pow4(Fonset1()))
                ,scalar(2)
            );
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Fonset3() const
{
    return max(scalar(1) - pow3(Rt()/scalar(2.5)),scalar(0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Fonset() const
{
    return max(Fonset2() - Fonset3(),scalar(0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::ReV() const
{
    return this->rho_*sqrt(2*magSqr(symm(fvc::grad(this->U()))))*sqr(this->y_)/this->mu();
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Rt() const
{
    return (this->rho_*this->k_)/(this->mu()*this->omega_);
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Flength() const
{
    tmp<volScalarField> Flength
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("Flength", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimless
        )
    );

    forAll(ReThetatTilda_, celli)
    {
    if (ReThetatTilda_[celli] < 400)
    {
        Flength.ref()[celli] = 398.189e-1 + (-119.270e-4)*ReThetatTilda_[celli] + (-132.567e-6)*ReThetatTilda_[celli];
    }
    else if (ReThetatTilda_[celli]< 596)
    {
        Flength.ref()[celli] = 263.404 + (-123.939e-2)*ReThetatTilda_[celli] + (-194.548e-5)*sqr(ReThetatTilda_[celli]) + (-101.695e-8)*pow3(ReThetatTilda_[celli]);
    }
    else if (ReThetatTilda_[celli] < 1200)
    {
        Flength.ref()[celli] = 0.5 - (ReThetatTilda_[celli] - 596.0)*3.0e-4;
    }
    else
    {
        Flength.ref()[celli] = 0.3188;
    }
    }   

    return Flength*(1 - Fsublayer()) + 40.0*Fsublayer();
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::ReThetac() const
{
   tmp<volScalarField> ReThetac
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("ReThetac", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimless
        )
    );
    forAll(ReThetatTilda_, celli)
    {
    if (ReThetatTilda_[celli] <= 1870)
    {
        ReThetac.ref()[celli] = ReThetatTilda_[celli] - (-396.035e-2 - 120.656e-4*ReThetatTilda_[celli] + 868.230e-6*sqr(ReThetatTilda_[celli]) + 696.506e-9*pow3(ReThetatTilda_[celli]) +174.105e-12*pow4(ReThetatTilda_[celli]));
    }
    else
    {
        ReThetac.ref()[celli] = ReThetatTilda_[celli] -(593.11 + (ReThetatTilda_[celli] -1870.0)*0.482);
    }
    }
    return ReThetac;
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Fsublayer() const
{
    return exp(-sqr(ReOmega()/200));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::ReOmega() const
{
    return this->rho_*this->omega_*sqr(this->y_)/this->mu();
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::gammaSep() const
{
    return min(
                2.0*max(
                        scalar(0),
                        ReV()/(3.235*ReThetac() - 1)*Freattach()
                    ),
                scalar(2)
            )*Fthetat();
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Freattach() const
{
    return exp(-pow4(Rt()/20));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::gammaEff() const
{
    return max(gamma_,gammaSep());
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Fturb() const
{
    return exp(-pow4(Rt()/4));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::PReThetatTilda() const
{
    // return tmp<volScalarField>
    // (
    //     new volScalarField
    //     (
    //         IOobject::groupName("PReThetatTilda", this->U_.group()),
    //         0.03*this->rho_/(500.0*this->mu()/(this->rho_*sqr(mag(this->U()))))*(ReThetat() - ReThetatTilda_)*(scalar(1.0) - Fthetat())
    //     )
    // );
    // return scalar(500)*(this->mu()/this->rho_)/(sqr(this->y_)*this->omega_); 
    return scalar(0.03)*this->rho_/(scalar(500.0)*this->mu()/(this->rho_*sqr(mag(this->U()))))*(ReThetat() - ReThetatTilda_)*(scalar(1.0) - Fthetat()) ;   
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Fthetat() const
{
    return min(
                max(
                    Fwake()*exp(-(this->y_/delta())),
                    1.0 - sqr((gamma_ - 1.0/50.0)/(1.0 - 1.0/50.0))
                    ),
            1.0
            );
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::ReThetat() const
{
    tmp<volScalarField> tmu = this->mu();
    volScalarField& mu = tmu.ref();

    tmp<volScalarField> ReThetat;
    // (
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             IOobject::groupName("ReThetat", this->U_.group()),
    //             this->runTime_.timeName(),
    //             this->mesh_
    //         ),
    //         this->mesh_,
    //         dimless
    //     )
    // );

    tmp<volScalarField> dUds = this->U() & (fvc::grad(this->U()) & this->U_)/magSqr(this->U());

    forAll(ReThetat.ref(), celli)
    {
        scalar lambda = 0.0;
        scalar thetat;
        int iter = 0;
        do
        {
            const scalar Tu
            (
                max(100*sqrt((2.0/3.0)*this->k_[celli])/mag(this->U_[celli]), scalar(0.027))
            );

            if (Tu <= 1.3)
            {
                thetat = (1173.51 - 589.428*Tu + 0.2196/sqr(Tu))*F(lambda, Tu, dUds.ref()[celli])*mu[celli]/(this->rho_[celli]*mag(this->U_[celli]));
            }
            else
            {
                thetat = 331.50*(pow(Tu - 0.5658, -0.671))*F(lambda, Tu, dUds.ref()[celli])*mu[celli]/(this->rho_[celli]*mag(this->U_[celli]));
            }

            lambda = this->rho_[celli]*sqr(thetat)*dUds.ref()[celli]/mu[celli];

            ++iter;
        } while (iter < 100);

        ReThetat.ref()[celli] = this->rho_[celli]*mag(this->U_[celli])*thetat/mu[celli];
        
    }

    return ReThetat;
}

template<class BasicTurbulenceModel>
scalar
transitionSST<BasicTurbulenceModel>::F(scalar lambda, scalar Tu, scalar dUds) const
{
    if (dUds <= scalar(0))
    {
        return 1.0 - (-21.986*lambda - 123.66*sqr(lambda) - 405.689*pow3(lambda))*exp(-pow(Tu/1.5,1.5));
    }
    else
    {
        return 1.0 + 0.275*(1 - exp(-35.0*lambda))*exp(-Tu/0.5);
    }
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::delta() const
{
    return 50*sqrt(2*magSqr(skew(fvc::grad(this->U()))))*this->y_/mag(this->U())
            *15.0/2.0
            *ReThetatTilda_*this->mu()/(this->rho_*mag(this->U()));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Fwake() const
{
    return exp(-sqr(ReOmega()/1.0e5));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::F1(const volScalarField& CDkOmega) const
{
    return max(kOmegaSST<BasicTurbulenceModel>::F1(CDkOmega),F3());
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::F3() const
{
    return exp(-pow(Ry()/scalar(120),8));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Ry() const
{
    return this->rho_*this->y_*sqrt(this->k_)/this->mu();
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::Pk(const volScalarField G) const
{   
    // volScalarField::Internal ge(gammaEff());

    // const tmp<volScalarField> ge(gammaEff());

    // const tmp<volScalarField::Internal> pk(kOmegaSST<BasicTurbulenceModel>::Pk(G));

    // tmp<volScalarField> pkRet;

    // forAll(pk,celli)
    // {
    //     pkRet[celli] = pk[celli]*ge[celli]
    // }

    // return pkRet;

    return kOmegaSST<BasicTurbulenceModel>::Pk(G)*gammaEff();
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
transitionSST<BasicTurbulenceModel>::epsilonByk
(
    const volScalarField& F1,
    const volScalarField& F2
) const
{
    return
        min(max(gammaEff(), scalar(0.1)), scalar(1))
       *kOmegaSST<BasicTurbulenceModel>::epsilonByk(F1, F2);
}

template<class BasicTurbulenceModel>
void
transitionSST<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    correctGammaReThetatTilda();
    kOmegaSST<BasicTurbulenceModel>::correct();
}

    // // Local references
    // const alphaField& alpha = this->alpha_;
    // const rhoField& rho = this->rho_;
    // const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    // const volVectorField& U = this->U_;
    // volScalarField& nut = this->nut_;
    // fv::options& fvOptions(fv::options::New(this->mesh_));
    // const volScalarField gammaEff_ = gammaEff();

    // BasicTurbulenceModel::correct();

    // volScalarField divU
    // (
    //     fvc::div(fvc::absolute(this->phi(), U))()()
    // );

    // tmp<volTensorField> tgradU = fvc::grad(U);
    // volScalarField S2(2*magSqr(symm(tgradU())));
    // volScalarField GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
    // volScalarField G(this->GName(), nut()*GbyNu);
    // tgradU.clear();

    // // Update omega and G at the wall
    // this->omega_.boundaryFieldRef().updateCoeffs();

    // volScalarField CDkOmega
    // (
    //     (2*this->alphaOmega2_)*(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    // );

    // volScalarField F1(this->F1(CDkOmega));
    // volScalarField F23(this->F23());

    // {
    //     volScalarField gamma(this->gamma(F1));
    //     volScalarField beta(this->beta(F1));

    //     // Turbulent frequency equation
    //     tmp<fvScalarMatrix> omegaEqn
    //     (
    //         fvm::ddt(alpha, rho, this->omega_)
    //       + fvm::div(alphaRhoPhi, this->omega_)
    //       - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_)
    //      ==
    //         alpha()*rho()*gamma
    //        *min
    //         (
    //             GbyNu,
    //             (this->c1_/this->a1_)*this->betaStar_*this->omega_()
    //            *max(this->a1_*this->omega_(), this->b1_*F23()*sqrt(S2()))
    //         )
    //       - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, this->omega_)
    //       - fvm::Sp(alpha()*rho()*beta*this->omega_(), this->omega_)
    //       - fvm::SuSp
    //         (
    //             alpha()*rho()*(F1() - scalar(1))*CDkOmega()/this->omega_(),
    //             this->omega_
    //         )
    //       + this->Qsas(S2(), gamma, beta)
    //       + this->omegaSource()
    //       + fvOptions(alpha, rho, this->omega_)
    //     );

    //     omegaEqn.ref().relax();
    //     fvOptions.constrain(omegaEqn.ref());
    //     omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
    //     solve(omegaEqn);
    //     fvOptions.correct(this->omega_);
    //     bound(this->omega_, this->omegaMin_);
    // }

    // // Turbulent kinetic energy equation
    // tmp<fvScalarMatrix> kEqn
    // (
    //     fvm::ddt(alpha, rho, this->k_)
    //   + fvm::div(alphaRhoPhi, this->k_)
    //   - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_)
    //  ==
    //     alpha()*rho()*this->Pk(G)*gammaEff_
    //   - ((2.0/3.0)*alpha()*rho()*divU, this->k_)
    //   - (alpha()*rho()*this->epsilonByk(F1, F23), this->k_*min(max(gammaEff_,scalar(0.1)),scalar(1.0)))
    //   + this->kSource()
    //   + fvOptions(alpha, rho, this->k_)
    // );

    // kEqn.ref().relax();
    // fvOptions.constrain(kEqn.ref());
    // solve(kEqn);
    // fvOptions.correct(this->k_);
    // bound(this->k_, this->kMin_);

    // this->correctNut(S2, F23);

template<class BasicTurbulenceModel>
void
transitionSST<BasicTurbulenceModel>::correctGammaReThetatTilda()
{

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    // const volVectorField& U = this->U_;
    // volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Transition momentum Reynolds number equation
    tmp<fvScalarMatrix> ReThetatTildaEqn
    (
        fvm::ddt(alpha, rho, ReThetatTilda_)
      + fvm::div(alphaRhoPhi, ReThetatTilda_)
      - fvm::laplacian(alpha*rho*DReThetatTildaEff(), ReThetatTilda_)
      ==
        alpha()*rho()*PReThetatTilda()
    );

    ReThetatTildaEqn.ref().relax();
    fvOptions.constrain(ReThetatTildaEqn.ref());
    solve(ReThetatTildaEqn);
    fvOptions.correct(ReThetatTilda_);
    // bound(ReThetatTilda_, this->ReThetatTildaMin_)

    // Intermittency equation
    tmp<fvScalarMatrix> gammaEqn
    (
        fvm::ddt(alpha, rho, gamma_)
      + fvm::div(alphaRhoPhi, gamma_)
      - fvm::laplacian(alpha*rho*DgammaEff(), gamma_)
      ==
        alpha()*rho()*Pgamma()
      + alpha()*rho()*scalar(0.06)*this->omega_*gamma_*Fturb()
      - fvm::Sp
      (
        alpha()*rho()*scalar(0.06)*this->omega_*gamma_*Fturb()*scalar(50),
        gamma_
      )
    );

    gammaEqn.ref().relax();
    fvOptions.constrain(gammaEqn.ref());
    solve(gammaEqn);
    fvOptions.correct(gamma_);
    // bound(gamma_, this->gammaMin_);

    // gammaEff_ = gammaEff();

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
