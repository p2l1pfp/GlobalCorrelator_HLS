library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package regionizer_data is 
    type particle is record
        pt : signed(15 downto 0);
        eta : signed(9 downto 0);
        phi : signed(9 downto 0);
        rest : std_logic_vector(27 downto 0);
    end record;

    subtype word64 is std_logic_vector(63 downto 0);

    function particle_to_w64(p : particle) return word64;
    function w64_to_particle(d : word64) return particle;
    function null_particle return particle;

    type particles is array(natural range <>) of particle;
    type w64s      is array(natural range <>) of word64;

    constant PHI_SHIFT : signed(9 downto 0) := to_signed(160, 10); -- 2*pi/9, size of a phi nonant, track finder sector or fiducial part of one PF region
    constant PHI_BORDER : signed(9 downto 0) := to_signed(58, 10); -- 0.25 (0.30 would be 69) 
    constant PHI_MARGIN_POS : signed(9 downto 0) := to_signed(+(160/2-58), 10);  -- half-width of fiducial MINUS border (half-size of gap are between sector N and sector N+2)
    constant PHI_MARGIN_NEG : signed(9 downto 0) := to_signed(-(160/2-58), 10);  -- same but with negative sign
    constant PHI_HALFWIDTH_POS : signed(9 downto 0) := to_signed(+(160/2+58), 10); -- half size of a full region (fiducial PLUS border)
    constant PHI_HALFWIDTH_NEG : signed(9 downto 0) := to_signed(-(160/2+58), 10);  

    constant PHI_CALOSHIFT    : signed(9 downto 0) := to_signed( 480,        10);  -- 2*pi/3, size of an HGCal sector
    constant PHI_CALOSHIFT1   : signed(9 downto 0) := to_signed( 320,        10);  -- 2*pi/3 - 2*pi/9, distance between center of hgcal sector 1 and pf region 1 = 2 * size of a phi nonant
    constant PHI_CALOEDGE_POS : signed(9 downto 0) := to_signed(+(480/2-58), 10);  -- +(half-size of calo sector)-border
    constant PHI_CALOEDGE_NEG : signed(9 downto 0) := to_signed(-(480/2-58), 10);  -- -(half-size of calo sector)+border


    constant PFII : natural := 4;
    constant NPFREGIONS : natural := 9;

    constant NTKSECTORS : natural := 9;
    constant NTKFIBERS : natural := 2;
    constant NTKFIFOS : natural := NTKFIBERS*3;
    constant NTKSORTED : natural := 24;
    constant NTKSTREAM : natural := (NTKSORTED+PFII-1)/PFII;

    constant NCALOSECTORS : natural := 3;
    constant NCALOFIBERS : natural := 4;
    constant NCALOFIFO0 : natural := NCALOFIBERS;
    constant NCALOFIFO12 : natural := 2*NCALOFIBERS;
    constant NCALOFIFOS : natural := NCALOFIFO0+2*NCALOFIFO12;
    constant NCALOSORTED : natural := 20;
    constant NCALOSTREAM : natural := (NCALOSORTED+PFII-1)/PFII;

end package;

package body regionizer_data is
    function particle_to_w64(p : particle) return word64 is
        variable ret : word64;
    begin
        ret( 9 downto  0) := std_logic_vector(p.eta);
        ret(19 downto 10) := std_logic_vector(p.phi);
        ret(47 downto 32) := std_logic_vector(p.pt);
        ret(31 downto 20) := (p.rest(11 downto  0));
        ret(63 downto 48) := (p.rest(27 downto 12));
        return ret;
    end particle_to_w64;

    function w64_to_particle(d : word64) return particle is
        variable ret : particle;
    begin
        ret.eta := signed(d( 9 downto  0));
        ret.phi := signed(d(19 downto 10));
        ret.pt  := signed(d(47 downto 32));
        ret.rest(11 downto  0) := (d(31 downto 20));
        ret.rest(27 downto 12) := (d(63 downto 48));
        return ret;
    end w64_to_particle;

    function null_particle return particle is
        variable ret : particle;
    begin
        ret.eta := to_signed(0, ret.eta'length);
        ret.phi := to_signed(0, ret.phi'length);
        ret.pt  := to_signed(0,  ret.pt'length);
        ret.rest := (others => '0');
        return ret;
    end null_particle;
end regionizer_data;


