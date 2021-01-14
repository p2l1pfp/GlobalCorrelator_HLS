library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
library unisim;
use unisim.vcomponents.all;

use work.regionizer_data.all;

entity fifo_merge3 is
    --generic(
    --    FIFO_INDEX : natural := 0
    --);
    port(
        ap_clk   : in std_logic;
        d1_in    : in particle;
        d2_in    : in particle;
        d3_in    : in particle;
        d1_valid : in std_logic;
        d2_valid : in std_logic;
        d3_valid : in std_logic;
        roll     : in  std_logic;
        d_out      : out particle;
        valid_out  : out std_logic;
        full1      : out std_logic;
        full2      : out std_logic;
        full3      : out std_logic;
        -- begin debug
        dbg_w64    : out std_logic_vector(63 downto 0);
        -- end debug
        roll_out   : out std_logic
    );
end fifo_merge3;

architecture Behavioral of fifo_merge3 is
    signal q2, q3 : particle;
    signal q2_valid, q3_valid : std_logic := '0';
    signal full2_i, full3_i : std_logic := '0';
begin

     logic: process(ap_clk) 
           variable load2, load3 : std_logic;
        begin
            if rising_edge(ap_clk) then
                if roll = '1' then
                    if d1_valid = '1' then
                        d_out <= d1_in;
                    elsif d2_valid = '1' then
                        d_out <= d2_in;
                    else
                        d_out <= d3_in;
                    end if;
                    valid_out <= d1_valid or d2_valid or d3_valid;
                    roll_out  <= '1';
                    q2       <= d2_in; 
                    q2_valid <= d1_valid and d2_valid;
                    q3       <= d3_in; 
                    q3_valid <= (d1_valid or d2_valid) and d3_valid;
                    full2_i <= d1_valid and d2_valid;
                    full3_i <= (d1_valid or d2_valid) and d3_valid;
                else
                    load2 := (d1_valid or q2_valid) and not full2_i;
                    load3 := (d1_valid or d2_valid or q2_valid or q3_valid) and not full3_i;
                    if d1_valid = '1' then
                        d_out <= d1_in; 
                    elsif q2_valid = '1' then
                        d_out <= q2;
                    elsif d2_valid = '1' then
                        d_out <= d2_in;
                    elsif q3_valid = '1' then
                        d_out <= q3;
                    else
                        d_out <= d3_in;
                    end if;
 
                    valid_out <= d1_valid or d2_valid or d3_valid or q2_valid or q3_valid;
                    roll_out  <= '0';
                    full1   <= '0';
                    full2_i <= d1_valid and (d2_valid or q2_valid);
                    full3_i <= (d1_valid or d2_valid or q2_valid) and (d3_valid or q3_valid);
                    if load2 = '1' then
                        q2 <= d2_in;
                        q2_valid <= d2_valid;
                    else
                        q2_valid <= d1_valid and q2_valid;
                    end if;
                    if load3 = '1' then
                        q3 <= d3_in;
                        q3_valid <= d3_valid;
                    else
                        q3_valid <= (d1_valid or d2_valid or q2_valid) and q3_valid;
                    end if;
                end if; 
            end if;
        end process;

        full2 <= full2_i;
        full3 <= full3_i;

        dbg_w64(15 downto 0) <= std_logic_vector(q2.pt);
        dbg_w64(16) <= q2_valid;
        dbg_w64(32 downto 17) <= std_logic_vector(q3.pt);
        dbg_w64(33) <= q3_valid;
        dbg_w64(63 downto 34) <= (others => '0');

end Behavioral;
