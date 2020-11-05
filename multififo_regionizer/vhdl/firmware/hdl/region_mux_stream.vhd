library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
library unisim;
use unisim.vcomponents.all;

use work.regionizer_data.all;

entity region_mux_stream is
    generic(
        NITEMS : natural;
        NREGIONS : natural;
        NSTREAM : natural;
        OUTII : natural
    );
    port(
        ap_clk  : in std_logic;
        roll    : in  std_logic;
        d_in     : in particles(NITEMS*NREGIONS-1 downto 0);
        valid_in : in std_logic_vector(NITEMS*NREGIONS-1 downto 0);
        d_out      : out particles(NSTREAM-1 downto 0);
        valid_out  : out std_logic_vector(NSTREAM-1 downto 0);
        roll_out   : out std_logic
    );
end region_mux_stream;

architecture Behavioral of region_mux_stream is
    signal regions : particles(NREGIONS*NITEMS-1 downto 0);
    signal valid   : std_logic_vector(NREGIONS*NITEMS-1 downto 0) := (others => '0');
    signal count   : integer range 0 to OUTII-1 := 0;
begin

     logic: process(ap_clk) 
           variable below : std_logic_vector(NITEMS-1 downto 0);
        begin
            if rising_edge(ap_clk) then
                if roll = '1' then
                    regions <= d_in;
                    valid   <= valid_in;
                    count    <= 0;
                else
                    if count < OUTII-1 then
                        regions(NITEMS-2 downto 0) <= regions(NITEMS-1 downto 1);
                        valid(NITEMS-2 downto 0)   <= valid(NITEMS-1 downto 1);
                        valid(NITEMS-1)            <= '0';
                        count <= count + 1;
                    else
                        for r in 0 to NREGIONS-2 loop
                            for i in 0 to NITEMS-1 loop
                                regions(r*NITEMS+i) <= regions((r+1)*NITEMS+i);
                                valid(  r*NITEMS+i) <= valid(  (r+1)*NITEMS+i);
                            end loop;
                        end loop;
                        for i in 0 to NITEMS-1 loop
                            valid((NREGIONS-1)*NITEMS+i) <= '0';
                        end loop;
                        count <= 0;
                    end if;
                end if;
                roll_out <= roll;
            end if;
        end process logic;

    gen_out: for i in 0 to NSTREAM-1 generate
        d_out(i) <= regions(i*OUTII);
        valid_out(i) <= valid(i*OUTII);
    end generate gen_out;


end Behavioral;
