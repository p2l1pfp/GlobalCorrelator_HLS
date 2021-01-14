library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.regionizer_data.all;

entity cascade_stream_sort is
    generic(
        NITEMS : natural := 24;
        NSTAGES : natural := 2
    );
    port(
        ap_clk  : in std_logic;
        roll    : in std_logic;
        d_in    : in particle;
        valid_in : in std_logic;
        shift_in : in std_logic := '0'; -- needed when cascading, leave to 0 otherwise
        d_out      : out particles(NITEMS-1 downto 0);
        valid_out  : out std_logic_vector(NITEMS-1 downto 0);
        roll_out   : out std_logic
    );
begin
    assert NSTAGES >= 2 report "BAD CHOICE OF NSTAGES = " & integer'image(NSTAGES);
end cascade_stream_sort;

architecture Behavioral of cascade_stream_sort is
    constant STAGE_ITEMS : natural := (NITEMS+NSTAGES-1)/NSTAGES;
    constant REST_ITEMS  : natural := NITEMS - STAGE_ITEMS;
    constant DELAY       : natural := NSTAGES-1;
    signal d_del     : particles(STAGE_ITEMS*DELAY-1 downto 0);
    signal valid_del : std_logic_vector(STAGE_ITEMS*DELAY-1 downto 0);
    signal roll_del  : std_logic := '0';
    signal d_carry : particle;
    signal valid_carry, shift_carry : std_logic := '0';
begin

    head : entity work.cascade_stream_sort_elem
                    generic map(NITEMS => STAGE_ITEMS)
                    port    map(ap_clk => ap_clk,
                                roll   => roll,
                                d_in   => d_in,
                                valid_in => valid_in,
                                shift_in => shift_in,
                                d_out     => d_del(STAGE_ITEMS*DELAY-1 downto STAGE_ITEMS*(DELAY-1)),
                                valid_out => valid_del(STAGE_ITEMS*DELAY-1 downto STAGE_ITEMS*(DELAY-1)),
                                roll_out  => open,
                                d_tail     => d_carry,
                                valid_tail => valid_carry,
                                shift_out  => shift_carry);
                                
    gen_tail_norec: if NSTAGES = 2 generate
        tail_norec : entity work.cascade_stream_sort_elem
                    generic map(NITEMS => REST_ITEMS)
                    port    map(ap_clk => ap_clk,
                                roll   => roll_del,
                                d_in   => d_carry,
                                valid_in => valid_carry,
                                shift_in => shift_carry,
                                d_out     => d_out(NITEMS-1 downto STAGE_ITEMS),
                                valid_out => valid_out(NITEMS-1 downto STAGE_ITEMS),
                                roll_out  => roll_out,
                                d_tail     => open,
                                valid_tail => open,
                                shift_out  => open);
    end generate gen_tail_norec;

    tail_rec: if NSTAGES > 2 generate
                tail_norec : entity work.cascade_stream_sort
                    generic map(NITEMS => REST_ITEMS, NSTAGES => NSTAGES-1)
                    port    map(ap_clk => ap_clk,
                                roll   => roll_del,
                                d_in   => d_carry,
                                valid_in => valid_carry,
                                shift_in => shift_carry,
                                d_out     => d_out(NITEMS-1 downto STAGE_ITEMS),
                                valid_out => valid_out(NITEMS-1 downto STAGE_ITEMS),
                                roll_out  => roll_out);
    end generate tail_rec;
         

    delay_logic: process(ap_clk) 
        begin
            if rising_edge(ap_clk) then
                roll_del <= roll;
                d_out(STAGE_ITEMS-1 downto 0) <= d_del(STAGE_ITEMS-1 downto 0);
                valid_out(STAGE_ITEMS-1 downto 0) <= valid_del(STAGE_ITEMS-1 downto 0);
                if DELAY > 1 then
                    d_del(STAGE_ITEMS*(DELAY-1)-1 downto 0) <= d_del(STAGE_ITEMS*DELAY-1 downto STAGE_ITEMS);
                    valid_del(STAGE_ITEMS*(DELAY-1)-1 downto 0) <= valid_del(STAGE_ITEMS*DELAY-1 downto STAGE_ITEMS);
                end if;
           end if;
        end process;

end Behavioral;
