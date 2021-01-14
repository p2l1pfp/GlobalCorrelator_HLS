library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
library unisim;
use unisim.vcomponents.all;

use work.regionizer_data.all;

entity stream_sort is
    generic(
        NITEMS : natural := 24
    );
    port(
        ap_clk  : in std_logic;
        roll    : in std_logic;
        d_in    : in particle;
        valid_in : in std_logic;
        d_out      : out particles(NITEMS-1 downto 0);
        valid_out  : out std_logic_vector(NITEMS-1 downto 0);
        roll_out   : out std_logic
    );
end stream_sort;

architecture Behavioral of stream_sort is
    signal sorted : particles(NITEMS-1 downto 0);
    signal valid  : std_logic_vector(NITEMS-1 downto 0) := (others => '0');
begin
     roll_out <= roll; -- the clock cycle a new event comes in here is also the clock cycle at which we're done with the old one

     logic: process(ap_clk) 
           variable below : std_logic_vector(NITEMS-1 downto 0);
        begin
            if rising_edge(ap_clk) then
                if roll = '1' then
                    sorted(0) <= d_in;
                    valid <= (0 => valid_in, others => '0'); 
                else
                    for i in NITEMS-1 downto 0 loop
                        if valid(i) = '0' or (valid_in = '1' and d_in.pt > sorted(i).pt) then
                            below(i) := '1';
                        else
                            below(i) := '0';
                        end if;
                    end loop;
                    for i in NITEMS-1 downto 1 loop
                        if below(i) = '1' and below(i-1) = '1' then
                            sorted(i) <= sorted(i-1);
                            valid(i)  <= valid(i-1);
                        elsif below(i) = '1' then
                            sorted(i) <= d_in;
                            valid(i)  <= valid_in;
                        else
                            sorted(i) <= sorted(i);
                            valid(i) <= valid(i);
                        end if;
                    end loop;
                    if below(0) = '1' then
                        sorted(0) <= d_in;
                        valid(0)  <= valid_in;
                    else
                        sorted(0) <= sorted(0);
                        valid(0) <= valid(0);
                    end if;
                end if; 
            end if;
        end process;

    d_out <= sorted;
    valid_out <= valid;

end Behavioral;
