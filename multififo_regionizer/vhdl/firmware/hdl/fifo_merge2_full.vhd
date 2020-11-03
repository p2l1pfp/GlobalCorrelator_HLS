library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
library unisim;
use unisim.vcomponents.all;

use work.regionizer_data.all;

entity fifo_merge2_full is
    --generic(
    --    FIFO_INDEX : natural := 0
    --);
    port(
        ap_clk   : in std_logic;
        d1_in    : in particle;
        d2_in    : in particle;
        d1_valid : in std_logic;
        d2_valid : in std_logic;
        roll     : in  std_logic;
        full     : in  std_logic;
        d_out      : out particle;
        valid_out  : out std_logic;
        full1      : out std_logic;
        full2      : out std_logic;
        -- begin debug
        dbg_w64    : out std_logic_vector(63 downto 0);
        -- end debug
        roll_out   : out std_logic
    );
end fifo_merge2_full;

architecture Behavioral of fifo_merge2_full is
    signal q1, q2, out_i : particle;
    signal q1_valid, q2_valid, valid_out_i, roll_out_i : std_logic := '0';
    signal full1_i, full2_i   : std_logic := '0';
begin

     logic: process(ap_clk) 
           variable load2 : std_logic;
        begin
            if rising_edge(ap_clk) then
                if roll = '1' then
                    if d1_valid = '1' then
                        out_i <= d1_in;
                    else
                        out_i <= d2_in;
                    end if;
                    valid_out_i <= d1_valid or d2_valid;
                    roll_out_i  <= '1';
                    q1       <= d1_in; 
                    q1_valid <= '0';
                    q2       <= d2_in; 
                    q2_valid <= d1_valid and d2_valid;
                    full1_i <= '0';
                    full2_i <= d1_valid and d2_valid;
                elsif roll_out_i = '0' and full = '1' and valid_out_i = '1' then
                    if full1_i = '0' then
                        q1       <= d1_in;
                        q1_valid <= d1_valid;
                        full1_i  <= d1_valid; -- don't signal full if a queue slot is available
                    end if;
                    if full2_i = '0' then
                        q2       <= d2_in;
                        q2_valid <= d2_valid;
                        full2_i  <= d2_valid;
                    end if;
                    roll_out_i <= '0';
                else
                    load2 := (d1_valid or q1_valid or q2_valid) and not (full2_i);
                    if q1_valid = '1' then
                        out_i <= q1; 
                    elsif d1_valid = '1' then
                        out_i <= d1_in; 
                    elsif q2_valid = '1' then
                        out_i <= q2;
                    else
                        out_i <= d2_in;
                    end if;
 
                    valid_out_i <= d1_valid or d2_valid or q1_valid or q2_valid;
                    roll_out_i <= '0';
                    full1_i <= '0';
                    full2_i <= (d1_valid or q1_valid) and (d2_valid or q2_valid);
                    if load2 = '1' then
                        q2 <= d2_in;
                        q2_valid <= d2_valid;
                    else
                        q2_valid <= ((d1_valid or q1_valid) and q2_valid);
                    end if;

                    q1_valid <= '0';
                end if; 
            end if;
        end process;

        full1 <= full1_i;
        full2 <= full2_i;
        d_out <= out_i;
        valid_out <= valid_out_i;
        roll_out <= roll_out_i;

        dbg_w64(15 downto 0) <= std_logic_vector(q2.pt);
        dbg_w64(16) <= q2_valid;
        dbg_w64(32 downto 17) <= std_logic_vector(q1.pt);
        dbg_w64(33) <= q1_valid;
        dbg_w64(63 downto 34) <= (others => '0');

end Behavioral;
