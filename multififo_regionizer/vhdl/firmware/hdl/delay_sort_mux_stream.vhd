library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.regionizer_data.all;

entity delay_sort_mux_stream is
    generic(
        NSORTED : natural;
        DELAY  : natural;
        NREGIONS : natural := NPFREGIONS;
        NSTREAM : natural;
        OUTII : natural := PFII240;
        SORT_NSTAGES : natural := 1
    );
    port(
        ap_clk : IN STD_LOGIC;
        d_in     : IN w64s(NREGIONS-1 downto 0);
        valid_in : IN STD_LOGIC_VECTOR(NREGIONS-1 downto 0);
        roll     : IN STD_LOGIC;
        d_out    : OUT w64s(NSTREAM-1 downto 0);
        roll_out : OUT STD_LOGIC
    );
end delay_sort_mux_stream;

architecture Behavioral of delay_sort_mux_stream is
    constant EFFDELAY : natural := DELAY-(SORT_NSTAGES-1);

    signal delayed:        w64s(EFFDELAY*NREGIONS-1 downto 0);
    signal delayed_valid:  std_logic_vector(EFFDELAY*NREGIONS-1 downto 0) := (others => '0');
    signal delayed_roll:   std_logic_vector(EFFDELAY-1 downto 0);

    signal sorted:        particles(NSORTED*NREGIONS-1 downto 0);
    signal sorted_valid:  std_logic_vector(NSORTED*NREGIONS-1 downto 0) := (others => '0');
    signal sorted_roll:   std_logic_vector(NREGIONS-1 downto 0) := (others => '0');

    signal mux :        particles(NSTREAM-1 downto 0);
    signal mux_valid :  std_logic_vector(NSTREAM-1 downto 0) := (others => '0');
    signal mux_roll :   std_logic := '0';

begin

    delayer : process(ap_clk)
        variable istart, iend: natural;
    begin
        if rising_edge(ap_clk) then
            istart := (EFFDELAY-1)*NREGIONS; iend := EFFDELAY*NREGIONS;
            delayed(      iend-1 downto istart) <= d_in(      NREGIONS-1 downto 0);
            delayed_valid(iend-1 downto istart) <= valid_in(NREGIONS-1 downto 0);
            delayed_roll(EFFDELAY-1)               <= roll;
            if EFFDELAY > 1 then
                delayed(      istart-1 downto 0) <= delayed(      iend-1 downto NREGIONS);
                delayed_valid(istart-1 downto 0) <= delayed_valid(iend-1 downto NREGIONS);
                delayed_roll(EFFDELAY-2 downto 0) <= delayed_roll(EFFDELAY-1 downto 1);
            end if;
        end if;
    end process delayer;

    gen_sorters: for isort in NREGIONS-1 downto 0 generate
        gen_sort_simple: if SORT_NSTAGES = 1 generate
            sorter : entity work.stream_sort
                        generic map(NITEMS => NSORTED)
                        port map(ap_clk => ap_clk,
                            d_in => w64_to_particle(delayed(isort)),
                            valid_in => delayed_valid(isort),
                            roll => delayed_roll(0),
                            d_out => sorted((isort+1)*NSORTED-1 downto isort*NSORTED),
                            valid_out => sorted_valid((isort+1)*NSORTED-1 downto isort*NSORTED),
                            roll_out => sorted_roll(isort)
                        );
            end generate gen_sort_simple;
        gen_sort_cascade: if SORT_NSTAGES > 1 generate
            sorter : entity work.cascade_stream_sort
                        generic map(NITEMS => NSORTED, NSTAGES => SORT_NSTAGES)
                        port map(ap_clk => ap_clk,
                            d_in => w64_to_particle(delayed(isort)),
                            valid_in => delayed_valid(isort),
                            roll => delayed_roll(0),
                            d_out => sorted((isort+1)*NSORTED-1 downto isort*NSORTED),
                            valid_out => sorted_valid((isort+1)*NSORTED-1 downto isort*NSORTED),
                            roll_out => sorted_roll(isort)
                        );
            end generate gen_sort_cascade;
        end generate gen_sorters;

    muxer: entity work.region_mux_stream
                    generic map(NREGIONS => NREGIONS, 
                                NITEMS   => NSORTED,
                                NSTREAM  => NSTREAM,
                                OUTII    => OUTII)
                    port map(ap_clk => ap_clk,
                        roll => sorted_roll(0),
                        d_in => sorted,
                        valid_in => sorted_valid,
                        d_out => mux,
                        valid_out => mux_valid,
                        roll_out => mux_roll);

     format: process(ap_clk)
        begin
            if rising_edge(ap_clk) then
                for i in 0 to NSTREAM-1 loop
                    if mux_valid(i) = '1' then
                        d_out(i) <= particle_to_w64(mux(i));
                    else
                        d_out(i) <= (others => '0');
                    end if;
                end loop;
                roll_out <= mux_roll;
            end if;
        end process format;

end Behavioral;
