{
  description = "EditR: a method to quantify base editing via Sanger sequencing";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/1559d3daa3ecc813a650b79375ea61b6741b8746";
  };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
      rpkgs = with pkgs.rPackages; [
        Biostrings
        pwalign
        sangerseqR
        gamlss
        gamlss_dist
        magrittr
        dplyr
        tidyr
        ggplot2
        cowplot
        gridExtra
        rmarkdown
        plotly
        yaml
        testthat
      ];
      editR = pkgs.rPackages.buildRPackage {
        name = "editR";
        src = self;
        propagatedBuildInputs = rpkgs;
      };
      rEnv = pkgs.rWrapper.override { packages = rpkgs ++ [ editR ]; };
      editrCli = pkgs.writeShellApplication {
        name = "editr";
        runtimeInputs = [ rEnv pkgs.pandoc ];
        text = ''
          Rscript ${self}/cli.R "$@"
        '';
      };
    in
    {
      packages.${system}.default = editR;

      devShells.${system}.default = pkgs.mkShell {
        buildInputs = [ rEnv pkgs.pandoc ];
      };

      apps.${system}.default = {
        type = "app";
        program = "${editrCli}/bin/editr";
      };
    };
}
