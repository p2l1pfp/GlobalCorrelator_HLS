library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.regionizer_data.all;

entity cascade_stream_sort_elem is
    generic(
        NITEMS : natural := 24
    );
    port(
        ap_clk  : in std_logic;
        roll    : in std_logic;
        d_in    : in particle;
        valid_in : in std_logic;
        shift_in : in std_logic;
        d_out      : out particles(NITEMS-1 downto 0);
        valid_out  : out std_logic_vector(NITEMS-1 downto 0);
        d_tail     : out particle;
        valid_tail : out std_logic;
        shift_out  : out std_logic;
        roll_out   : out std_logic
    );
end cascade_stream_sort_elem;


architecture Behavioral of cascade_stream_sort_elem is
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
                        if shift_in = '1' or (below(i) = '1' and below(i-1) = '1') then
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
                    if  shift_in = '1' or below(0) = '1' then
                        sorted(0) <= d_in;
                        valid(0)  <= valid_in;
                    else
                        sorted(0) <= sorted(0);
                        valid(0) <= valid(0);
                    end if;
                    if shift_in = '1' or below(NITEMS-1) = '1' then
                        d_tail     <= sorted(NITEMS-1);
                        valid_tail <=  valid(NITEMS-1);
                        shift_out  <=  valid(NITEMS-1);
                    else
                        d_tail     <= d_in;
                        valid_tail <= valid_in;
                        shift_out  <= '0';
                    end if;
                end if; 
            end if;
        end process;

    d_out <= sorted;
    valid_out <= valid;

end Behavioral;

