library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
library unisim;
use unisim.vcomponents.all;

use work.regionizer_data.all;

entity fifo_merge2 is
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
        --out_full : in  std_logic;
        d_out      : out particle;
        valid_out  : out std_logic;
        full1      : out std_logic;
        full2      : out std_logic;
        -- begin debug
        dbg_w64    : out std_logic_vector(63 downto 0);
        -- end debug
        roll_out   : out std_logic
    );
end fifo_merge2;

architecture Behavioral of fifo_merge2 is
    signal queue : particle;
    signal queue_valid : std_logic := '0';
    signal full2_i     : std_logic := '0';
begin

     logic: process(ap_clk) 
           variable load2 : std_logic;
        begin
            if rising_edge(ap_clk) then
                if roll = '1' then
                    if d1_valid = '1' then
                        d_out <= d1_in;
                    else
                        d_out <= d2_in;
                    end if;
                    valid_out <= d1_valid or d2_valid;
                    roll_out  <= '1';
                    queue      <= d2_in; 
                    queue_valid <= d1_valid and d2_valid;
                    full1 <= '0';
                    full2_i <= d1_valid and d2_valid;
                else
                    load2 := (d1_valid or queue_valid) and not (full2_i);
                    if d1_valid = '1' then
                        d_out <= d1_in; 
                    elsif queue_valid = '1' then
                        d_out <= queue;
                    else
                        d_out <= d2_in;
                    end if;
 
                    valid_out <= d1_valid or d2_valid or queue_valid;
                    roll_out <= '0';
                    full1 <= '0';
                    full2_i <= d1_valid and (d2_valid or queue_valid);
                    if load2 = '1' then
                        queue <= d2_in;
                        queue_valid <= d2_valid;
                    else
                        queue_valid <= (d1_valid and queue_valid); -- queue_valid when (d1_valid and queue_valid) else '0';
                    end if;
                end if; 
            end if;
        end process;

        full2 <= full2_i;

        dbg_w64(15 downto 0) <= std_logic_vector(queue.pt);
        dbg_w64(16) <= queue_valid;
        dbg_w64(63 downto 17) <= (others => '0');

end Behavioral;
